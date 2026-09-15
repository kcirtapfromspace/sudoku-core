use super::*;

#[test]
fn test_empty() {
    let set = BitSet::empty();
    assert!(set.is_empty());
    assert_eq!(set.count(), 0);
}

#[test]
fn test_all_9() {
    let set = BitSet::all_9();
    assert_eq!(set.count(), 9);
    for i in 1..=9 {
        assert!(set.contains(i));
    }
    assert!(!set.contains(0));
    assert!(!set.contains(10));
}

#[test]
fn test_insert_remove() {
    let mut set = BitSet::empty();
    set.insert(5);
    assert!(set.contains(5));
    assert_eq!(set.count(), 1);

    set.insert(3);
    assert!(set.contains(3));
    assert_eq!(set.count(), 2);

    set.remove(5);
    assert!(!set.contains(5));
    assert_eq!(set.count(), 1);
}

#[test]
fn test_single_value() {
    let mut set = BitSet::empty();
    set.insert(7);
    assert_eq!(set.single_value(), Some(7));

    set.insert(3);
    assert_eq!(set.single_value(), None);
}

#[test]
fn test_iter() {
    let set = BitSet::from_slice(&[1, 3, 5, 9]);
    let values: Vec<u8> = set.iter().collect();
    assert_eq!(values, vec![1, 3, 5, 9]);
}

#[test]
fn test_operations() {
    let a = BitSet::from_slice(&[1, 2, 3]);
    let b = BitSet::from_slice(&[2, 3, 4]);

    let union = a.union(&b);
    assert_eq!(union.to_vec(), vec![1, 2, 3, 4]);

    let intersection = a.intersection(&b);
    assert_eq!(intersection.to_vec(), vec![2, 3]);

    let difference = a.difference(&b);
    assert_eq!(difference.to_vec(), vec![1]);
}
