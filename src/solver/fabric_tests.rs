use super::*;

#[test]
fn test_cell_index_roundtrip() {
    for row in 0..9 {
        for col in 0..9 {
            let idx = cell_index(row, col);
            assert_eq!(cell_pos(idx), (row, col));
        }
    }
}

#[test]
fn test_sector_cells() {
    // Row 0
    let row0 = sector_cells(0);
    assert_eq!(row0, [0, 1, 2, 3, 4, 5, 6, 7, 8]);

    // Col 0
    let col0 = sector_cells(9);
    assert_eq!(col0, [0, 9, 18, 27, 36, 45, 54, 63, 72]);

    // Box 0
    let box0 = sector_cells(18);
    assert_eq!(box0, [0, 1, 2, 9, 10, 11, 18, 19, 20]);
}

#[test]
fn test_peers() {
    let peers = compute_peers(0); // cell (0,0)
    assert_eq!(peers.len(), 20);
    // Should include (0,1)..(0,8), (1,0)..(8,0), and box peers not in row/col
    assert!(peers.contains(&1)); // (0,1)
    assert!(peers.contains(&9)); // (1,0)
    assert!(peers.contains(&10)); // (1,1) - box peer
}

#[test]
fn test_fabric_from_grid() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();
    let fab = CandidateFabric::from_grid(&grid);

    // Cell (0,0) = 5, should be a given
    assert_eq!(fab.values[0], Some(5));
    assert!(fab.is_given[0]);

    // Cell (0,2) is empty, should have candidates
    let idx = cell_index(0, 2);
    assert!(fab.values[idx].is_none());
    assert!(!fab.cell_cands[idx].is_empty());

    // Verify candidate doesn't include 5 (same row as (0,0))
    assert!(!fab.cell_cands[idx].contains(5));
}

#[test]
fn test_sees() {
    let puzzle =
        "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
    let grid = Grid::from_string(puzzle).unwrap();
    let fab = CandidateFabric::from_grid(&grid);

    // Same row
    assert!(fab.sees(0, 5));
    // Same col
    assert!(fab.sees(0, 9));
    // Same box
    assert!(fab.sees(0, 10));
    // Not seeing each other
    assert!(!fab.sees(0, 40)); // (0,0) and (4,4)
}
