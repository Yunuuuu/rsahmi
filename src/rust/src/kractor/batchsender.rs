use crossbeam_channel::{SendError, Sender};

const DEFEALT_BATCH_SIZE: usize = 20;

#[derive(Clone)]
pub struct BatchSender<T> {
    msg_vec: Vec<T>,
    tx: Sender<Vec<T>>,
    capacity: usize,
}

impl<T> BatchSender<T> {
    #[allow(dead_code)]
    pub fn new(sender: Sender<Vec<T>>) -> Self {
        Self::with_capacity(DEFEALT_BATCH_SIZE, sender)
    }

    pub fn with_capacity(capacity: usize, sender: Sender<Vec<T>>) -> Self {
        Self {
            msg_vec: Vec::with_capacity(capacity),
            tx: sender,
            capacity,
        }
    }

    pub fn send(&mut self, msg: T) -> Result<(), SendError<Vec<T>>> {
        if self.msg_vec.len() >= self.msg_vec.capacity() {
            let mut pack = Vec::with_capacity(self.capacity);
            std::mem::swap(&mut self.msg_vec, &mut pack);
            self.tx.send(pack)?
        }
        self.msg_vec.push(msg);
        Ok(())
    }

    pub fn flush(&mut self) -> Result<(), SendError<Vec<T>>> {
        if !self.msg_vec.is_empty() {
            self.tx.send(std::mem::take(&mut self.msg_vec))?;
        }
        Ok(())
    }
}

impl<T> Drop for BatchSender<T> {
    fn drop(&mut self) {
        if !self.msg_vec.is_empty() {
            // Just omit the Error message
            let _ = self.tx.send(std::mem::take(&mut self.msg_vec));
        }
    }
}

#[cfg(test)]
mod tests {
    use crossbeam_channel::unbounded;

    use super::*;

    #[test]
    fn test_send_less_than_capacity() {
        let (tx, rx) = unbounded();
        let mut batcher = BatchSender::with_capacity(3, tx);
        batcher.send(1).unwrap();
        batcher.send(2).unwrap();
        assert!(rx.try_recv().is_err()); // Nothing sent yet
    }

    #[test]
    fn test_send_more_capacity() {
        let (tx, rx) = unbounded();
        let mut batcher = BatchSender::with_capacity(3, tx);
        batcher.send(1).unwrap();
        batcher.send(2).unwrap();
        batcher.send(3).unwrap();
        batcher.send(4).unwrap(); // should trigger send
        let batch = rx.try_recv().unwrap();
        assert_eq!(batch, vec![1, 2, 3]);
    }

    #[test]
    fn test_send_multiple_batches() {
        let (tx, rx) = unbounded();
        let mut batcher = BatchSender::with_capacity(2, tx);
        for i in 0 .. 5 {
            batcher.send(i).unwrap();
        }
        let batch1 = rx.try_recv().unwrap();
        let batch2 = rx.try_recv().unwrap();
        assert_eq!(batch1, vec![0, 1]);
        assert_eq!(batch2, vec![2, 3]);
        assert!(rx.try_recv().is_err()); // 4 is still buffered
    }

    #[test]
    fn test_explicit_flush() {
        let (tx, rx) = unbounded();
        let mut batcher = BatchSender::with_capacity(5, tx);
        batcher.send(42).unwrap();
        batcher.flush().unwrap();
        let batch = rx.try_recv().unwrap();
        assert_eq!(batch, vec![42]);
    }

    #[test]
    fn test_drop_flushes() {
        let (tx, rx) = unbounded();
        {
            let mut batcher = BatchSender::with_capacity(10, tx);
            batcher.send(99).unwrap();
        } // dropped
        let batch = rx.recv().unwrap();
        assert_eq!(batch, vec![99]);
    }
}
