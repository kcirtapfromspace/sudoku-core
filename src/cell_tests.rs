use super::*;

#[test]
fn test_new_empty() {
    let cell = Cell::new_empty();
    assert!(cell.is_empty());
    assert!(!cell.is_given());
    assert_eq!(cell.candidate_count(), 9);
}

#[test]
fn test_new_given() {
    let cell = Cell::new_given(5);
    assert!(cell.is_filled());
    assert!(cell.is_given());
    assert_eq!(cell.value(), Some(5));
    assert_eq!(cell.candidate_count(), 0);
}

#[test]
fn test_clear_given() {
    let mut cell = Cell::new_given(5);
    cell.clear();
    // Given cells should not be clearable
    assert!(cell.is_filled());
    assert_eq!(cell.value(), Some(5));
}

#[test]
fn test_clear_non_given() {
    let mut cell = Cell::new_filled(5);
    cell.clear();
    assert!(cell.is_empty());
    assert_eq!(cell.candidate_count(), 9);
}

#[test]
fn test_candidates() {
    let mut cell = Cell::new_empty();
    cell.remove_candidate(3);
    cell.remove_candidate(7);
    assert!(!cell.has_candidate(3));
    assert!(!cell.has_candidate(7));
    assert!(cell.has_candidate(5));
    assert_eq!(cell.candidate_count(), 7);

    cell.toggle_candidate(3);
    assert!(cell.has_candidate(3));
}
