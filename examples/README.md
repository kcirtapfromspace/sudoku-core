# Sudoku Engine Examples

## Basic Usage

Run the basic example to see puzzle generation and solving:

```bash
cargo run --example basic
```

## TUI Application

Run the interactive terminal interface:

```bash
cargo run --bin sudoku
```

### Controls

| Key | Action |
|-----|--------|
| Arrow keys / hjkl | Navigate cells |
| w/a/s/d | Jump to adjacent 3x3 box |
| 1-9 | Enter number |
| Shift + 1-9 | Toggle candidate |
| 0 / Delete / Backspace | Clear cell |
| c | Toggle candidate mode |
| f | Fill candidates (current cell) |
| Shift+F | Fill candidates (all cells) |
| x | Clear notes (current cell) |
| Shift+X | Clear notes (all cells) |
| g | Toggle ghost hints |
| v | Toggle valid cell highlighting |
| u | Undo |
| Ctrl+r | Redo |
| ? | Show hint (press again for proof detail) |
| ! | Apply hint |
| n | New game menu |
| p | Pause |
| Shift+S | Stats screen |
| q | Quit |

## Building for Mobile

### iOS (Swift)

```bash
cd crates/sudoku-ffi
cargo build --release --target aarch64-apple-ios
cargo run --bin uniffi-bindgen generate --library ../target/release/libsudoku_ffi.dylib --language swift --out-dir ./generated
```

### Android (Kotlin)

```bash
cd crates/sudoku-ffi
cargo build --release --target aarch64-linux-android
cargo run --bin uniffi-bindgen generate --library ../target/release/libsudoku_ffi.so --language kotlin --out-dir ./generated
```
