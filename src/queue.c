#include "queue.h"

typedef struct queueNodeTag
{
  queueElementT element;
  struct queueNodeTag *next;
} queueNodeT;

typedef struct queueCDT
{
  queueNodeT *front, *rear;
  int length;
} queueCDT;

queueADT QueueCreat(void)
{
  queueADT queue;
  
  queue = (queueADT) malloc(sizeof(queue));

  if(queue == NULL)
    {
      fprintf(stderr, "Insufficient momory for new queue.\n");
      exit(1); /* Exit program, returning error code.*/
    }
  queue->front = queue->rear = NULL;
  
  queue->length = 0;

  return queue;
}

void QueueDestroy(queueADT queue)
{
  /* 
   * First remove each element from the queue (each elemet
   * is in a dynamically-allocated node.)
   */
  while(!QueueIsEmpty(queue))
    QueueDelete(queue);

  /*
   * Reset the front and rear just in case someone 
   * tries to use them after the CDT is freed.
   */
  queue->front = queue->rear = NULL;
  queue->length = 0;

  /*
   * Now free the structure that holds information
   * about the queue.
   */
  free(queue);
}

void QueueEnter(queueADT queue, queueElementT element)
{
  queueNodeT *newNodeP;
  /* Allocate space for a node in the linked list. */

  newNodeP = (queueNodeT *)malloc(sizeof(queueNodeT));
  
  if(newNodeP == NULL)
    {
      fprintf(stderr,"Insufficient memory for new queue element.\n");
      exit(0); /*Exit program, returnning error code.*/
    }
  
  /*Place information in the node.*/
  newNodeP->element = element;
  newNodeP->next = NULL;
  
  /* 
   * Link the element into the right place in 
   * the linked list.
   */
  
  if(queue->front == NULL)/* queue is empty */
    {
      queue->front = queue->rear = newNodeP;
      queue->length = 1;
    }
  else
    {
      queue->rear->next = newNodeP;
      queue->rear = newNodeP;
      queue->length += 1;
    }
}

queueElementT QueueDelete(queueADT queue)
{
  if(queue->front == NULL)
    { 
      fprintf(stderr,"Try deleting a node from an empty queue.\n");
      exit(0); /*Exit program, returnning error code.*/
    }

  queueElementT tmp = queue->front->element;
  queueNodeT *newfront = queue->front->next;
  
  free(queue->front);
  queue->front = newfront;
  queue->length -= 1;
  
  return tmp;
}

int QueueIsEmpty(queueADT queue)
{
  if(queue->front ==NULL)
    return 1;
  else 
    return 0;
}

int QueueLength(queueADT queue)
{
  return queue->length;
}
