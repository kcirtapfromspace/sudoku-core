#!/usr/bin/env node
// Parse every authored Mermaid block with the pinned documentation tools.
import { readFile } from 'node:fs/promises';
import { execFileSync } from 'node:child_process';
import { resolve, dirname } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const dependencies = resolve(root, process.argv[2] ?? 'target/docs-tools/node_modules');
const { JSDOM } = await import(pathToFileURL(resolve(dependencies, 'jsdom/lib/api.js')));
const dom = new JSDOM('<!doctype html><html><body></body></html>');
globalThis.window = dom.window;
globalThis.document = dom.window.document;
const { default: mermaid } = await import(pathToFileURL(resolve(dependencies, 'mermaid/dist/mermaid.esm.mjs')));
mermaid.initialize({ startOnLoad: false, securityLevel: 'strict' });

const files = [...new Set(execFileSync('git', [
  'ls-files', '--cached', '--others', '--exclude-standard', '--', '*.md', '*.mdx',
], { cwd: root, encoding: 'utf8' }).trim().split('\n').filter(Boolean))];
let count = 0;
let failed = false;
for (const file of files) {
  const markdown = await readFile(resolve(root, file), 'utf8');
  const blocks = [];
  let fence;
  const lines = markdown.split('\n');
  for (let i = 0; i < lines.length; i += 1) {
    const marker = lines[i].match(/^ {0,3}(`{3,}|~{3,})(.*)$/);
    if (!fence && marker) {
      fence = { marker: marker[1], language: marker[2].trim(), line: i + 1, body: [] };
    } else if (fence && marker && marker[1][0] === fence.marker[0]
      && marker[1].length >= fence.marker.length && !marker[2].trim()) {
      if (fence.language === 'mermaid') blocks.push(fence);
      fence = undefined;
    } else if (fence) {
      fence.body.push(lines[i]);
    }
  }
  if (fence) {
    failed = true;
    console.error(`${file}:${fence.line}: Unclosed code fence.`);
  }
  for (const block of blocks) {
    count += 1;
    try {
      await mermaid.parse(block.body.join('\n'));
    } catch (error) {
      failed = true;
      console.error(`${file}:${block.line}: ${error.message}`);
    }
  }
}
dom.window.close();
console.log(`Checked ${count} Mermaid diagrams in ${files.length} documents.`);
process.exitCode = failed ? 1 : 0;
