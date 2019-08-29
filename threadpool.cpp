#include <threadpool.h>

void threadpool::start(void) {
  for ( unsigned i=0; i<this->poolSize; i++ )
    pthread_create(&this->tid[i],NULL, threadpool::thread, this);
  return;
}

void threadpool::stop(void) {

  this->shutting_down = true;
  pthread_cond_broadcast(&this->queue_not_empty);

  for ( unsigned i=0; i<this->poolSize; i++ ) {
    if (this->verbose ) cout << "Calling cancel for " << this->tid[i] << "\n";
    pthread_cancel(this->tid[i]);
  }
  return;
}

void threadpool::set_thread( void fptr(void *) ) {
  this->funcPtr = fptr;
  return;
}

threadpool threadpool::operator++(int) {
  threadpool temp = *this;
  this->increase_pool();
  return temp;
}
threadpool threadpool::operator--(int) {
  threadpool temp = *this;
  this->decrease_pool();
  return temp;
}
threadpool& threadpool::operator++(void) {
  this->increase_pool();
  return *this;
}
threadpool& threadpool::operator--(void) {
  this->decrease_pool();
  return *this;
}

threadpool::threadpool( void fptr(void *), unsigned n, bool vb ) {

  tid = new pthread_t[n];

  queue_not_empty = queue_is_empty = PTHREAD_COND_INITIALIZER;
  print_lock = write_lock = read_lock = PTHREAD_MUTEX_INITIALIZER;

  queue_empty = true;
  verbose = vb;

  funcPtr = fptr;

  thread_id = 0;

  for ( unsigned i=0; i<n; i++ ) {
    pthread_create(&tid[i],NULL, threadpool::thread, this);
    usleep(10);
  }
  poolSize = n;
  shutting_down = false;

  return;
}

threadpool::~threadpool(void) {
  if (this->verbose ) cout << this << " shutting down\n";

  this->shutting_down = true;

  // Wake up the echoes... wait, no, wake up the threads
  // and let them know it's time to bail.
  pthread_cond_broadcast(&this->queue_not_empty);

  for ( unsigned i=0; i<this->poolSize; i++ ) {
    if (this->verbose ) cout << "Calling cancel for " << tid[i] << "\n";
    pthread_cancel(tid[i]);
    /*if (this->verbose ) cout << "Waiting for join\n";
    pthread_join(tid[i],retval);
    if (this->verbose ) cout << tid[i] << " complete\n";*/
  }

  return;
}

// Lock down the queue (dis-allow dequeuing) while enqueuing data
void threadpool::queue_lock( void ) {
  pthread_mutex_lock(&this->read_lock);
  return;
}

void threadpool::queue_unlock( void ) {
  pthread_mutex_unlock(&this->read_lock);
  return;
}

// Put something (a pointer to your data) onto the queue for processing
void threadpool::enqueue( void *p ) {

  pthread_mutex_lock(&this->write_lock);

  if ( this->verbose ) cout << "\nEQ: enqueuing\n";

  this->queue.push_back(p);
  this->queue_empty = false;

  if ( this->verbose ) cout << "\nEQ: Done queue has " << this->queue.size() << " elements\n";

  if ( this->verbose ) cout << "\nEQ: releasing lock\n";
  pthread_mutex_unlock(&this->write_lock);
  pthread_cond_broadcast(&this->queue_not_empty);

  return;
}

// Pull the next item off the queue and return it to the processing thread
void *threadpool::dequeue(pthread_t id) {
  void *p = NULL;
  void *last = NULL;

  if ( this->shutting_down ) {
    if ( this->verbose ) cout << id << " returning null on shutdown\n";
    return p;
  }
  
  if (this->verbose ) cout << "\n" << id << ": Waiting for read lock\n";
  pthread_mutex_lock(&this->read_lock);

  if ( this->queue_empty ) {
    if (this->verbose ) cout << "\n" << id << ": Waiting for queue fill\n";
    pthread_cond_wait(&this->queue_not_empty, &read_lock);
  }

  if (this->verbose ) cout << "\n" << id << ": Got lock... dequeuing\n";

  if ( this->queue.size() ) {

    p = (void *)0x0;

    p = this->queue.front();
    this->queue.pop_front();

    while ( !p || p == last ) {
      p = this->queue.front();
      this->queue[0] = (void *)0x0;
      this->queue.pop_front();
    }

    last = p;

    if ( this->verbose )
      cout << id << " unshifted releasing queue lock\n";

  } else {
    if (this->verbose ) cout << "\n" << id << ": Empty queue.\n";
    this->queue_empty = true;
  }

  if ( this->queue_empty )
    pthread_cond_broadcast(&this->queue_is_empty);

  pthread_mutex_unlock(&this->read_lock);

  if (this->verbose ) cout << "\n" << id << ": Done.\n";

  return p;

}

// Block until the queue is empty
void threadpool::wait_until_empty(void) {
  pthread_mutex_lock(&this->wait_lock);
  if ( !this->queue_empty ) {
    pthread_cond_wait(&this->queue_is_empty, &this->wait_lock);
  }
  pthread_mutex_unlock(&this->wait_lock);
  return;
}

// Utility functions ///////////////////////////////////////////////////////////////////////

unsigned int threadpool::get_pool_size(void) { return this->poolSize; }

unsigned int threadpool::get_queue_size(void) { return this->queue.size(); }

void threadpool::dump_queue(unsigned int items) {
  if ( !items )
    items = this->queue.size();
  for ( uint i=0; i<items; i++ )
    cout << this->queue[i] << endl;
  return;
}

void threadpool::increase_pool(unsigned int newThreads) {

  for ( unsigned int i=0; i<newThreads; i++ ) {
    pthread_create(&this->tid[this->poolSize],NULL, threadpool::thread, this);
    this->poolSize++;
  }
  return;
}

void threadpool::decrease_pool(unsigned int loose_threads) {
  unsigned int new_size = (this->poolSize-loose_threads > 0) ? this->poolSize-loose_threads : 1;

  if ( new_size < this->poolSize ) {
    for ( unsigned int i=this->poolSize-1; i>=new_size; i-- ) {
      pthread_cancel(this->tid[i]);
      this->poolSize--;
    }
  }
  return;
}
