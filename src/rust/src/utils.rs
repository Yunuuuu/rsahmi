pub unsafe fn split_chunk_inclusive<T>(
    chunk: &mut Vec<T>,
    pos: usize,
) -> Vec<T> {
    let len = chunk.len();
    let cap = chunk.capacity();
    let ptr = chunk.as_mut_ptr();
    let head_len = pos + 1;
    let tail_len = len - head_len;
    let tail_cap = cap - head_len;
    let tail_ptr = ptr.add(head_len);

    // Take ownership of chunk without dropping it on scope exit
    let _original = std::mem::ManuallyDrop::new(std::mem::take(chunk));

    // Create head Vec from start to pos (inclusive)
    let head = Vec::from_raw_parts(ptr, head_len, head_len);

    // Create tail Vec from pos+1 to end
    let tail = Vec::from_raw_parts(tail_ptr, tail_len, tail_cap);

    // Now assign tail back to chunk
    *chunk = tail;

    head
}

#[test]
fn test_split_chunk_inclusive() {
    let mut data = vec![1, 2, 3, 4, 5];
    let head = unsafe { split_chunk_inclusive(&mut data, 2) };

    // `head` should contain elements 1, 2, 3
    assert_eq!(head, vec![1, 2, 3]);

    // `data` should now contain the remaining elements 4, 5
    assert_eq!(data, vec![4, 5]);
}
