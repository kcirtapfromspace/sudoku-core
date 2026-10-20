// Run after generating the offline example. Uses the documentation's pinned JSDOM.
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import {pathToFileURL} from 'node:url';

const html = path.resolve(process.argv[2] ?? 'target/grid-viewer.html');
const dependencies = path.resolve(process.argv[3] ?? 'target/docs-tools/node_modules');
const {JSDOM, VirtualConsole} = await import(pathToFileURL(path.join(dependencies, 'jsdom/lib/api.js')));
const errors = [];
const console = new VirtualConsole();
console.on('jsdomError', error => errors.push(String(error)));
const dom = new JSDOM(fs.readFileSync(html, 'utf8'), {runScripts: 'dangerously', virtualConsole: console});
const {document, Event} = dom.window;
const get = id => document.getElementById(id);
const entries = JSON.parse(get('scene-data').textContent);
assert(entries.length > 0);
assert.equal(get('examples').children.length, entries.length);
const certified = entries.filter(e => e.proof_bundle);
assert.equal(certified.length, 48);
assert.equal(new Set(certified.map(e => e.scene.technique)).size, 48);
assert(certified.every(e => e.scene.evidence.status === 'checked_deduction'));
let frames = 0;
for (let i = 0; i < entries.length; i++) {
  get('examples').children[i].click();
  const {scene} = entries[i];
  assert.equal(get('example-title').textContent, entries[i].title);
  assert(get('status').classList.contains(scene.evidence.status));
  assert.equal(get('frame-select').options.length, scene.frames.length);
  assert.equal(get('export-proof').disabled, !entries[i].proof_bundle);
  for (let frame = 0; frame < scene.frames.length; frame++) {
    get('frame-select').value = frame;
    get('frame-select').dispatchEvent(new Event('input'));
    assert.equal(get('frame-count').textContent, `${frame + 1} / ${scene.frames.length}`);
    assert.equal(get('frame-title').textContent, scene.frames[frame].title);
    assert.equal(get('grid').querySelectorAll('[role=button]').length, 81);
    assert.equal(get('grid').querySelectorAll('line').length, 20 + scene.frames[frame].links.length);
    const cell = get('grid').querySelector('[role=button]');
    cell.dispatchEvent(new dom.window.MouseEvent('click'));
    assert.equal(get('cell-title').textContent, 'r1c1');
    frames++;
  }
  if (entries[i].original_puzzle.length === 81) {
    get('original').checked = true;
    get('original').dispatchEvent(new Event('change'));
    assert(get('board-mode').textContent.includes('Original givens'));
    assert.equal(get('grid').querySelectorAll('line').length, 20);
    get('original').checked = false;
    get('original').dispatchEvent(new Event('change'));
  } else {
    assert(get('original').disabled);
  }
  if (scene.frames.length > 1) {
    get('previous').click();
    assert.equal(get('frame-select').value, String(scene.frames.length - 2));
    get('next').click();
    assert.equal(get('frame-select').value, String(scene.frames.length - 1));
  }
}
get('search').value = 'no_such_example_315';
get('search').dispatchEvent(new Event('input'));
assert.equal(get('examples').children.length, 0);
get('search').value = '';
get('search').dispatchEvent(new Event('input'));
assert.equal(get('examples').children.length, entries.length);
assert.deepEqual(errors, []);
dom.window.close();
process.stdout.write(`Checked ${entries.length} examples and ${frames} grid frames.\n`);
