/*
 * Copyright 1993-2002 The MathWorks, Inc.
 * $Revision: 336 $  $Date: 2002/03/15 15:28:59 $
 */

#ifndef QUEUE_H
#define QUEUE_H
#include <stdio.h>
#include <stdlib.h>

typedef int queueElementT;

typedef struct queueCDT *queueADT;

queueADT QueueCreat();
void QueueDestroy(queueADT queue);
int QueueLength(queueADT queue);
void QueueEnter(queueADT queue, queueElementT item);
queueElementT QueueDelete(queueADT queue);
int QueueIsEmpty(queueADT queue);

#endif /* QUEUE_H */
