 // Header file for high efficiency banded linear equation system solver (BLESS)
 // Author: William Read, Universidad Nacional Pedro Henríquez Ureña, Santo Domingo.
 // May be used, copied, reproduced or modified, giving credit to the author.
 //
 
#ifndef BLESS_H
#define BLESS_H
#endif

enum task { trian, reslv
};

// Declare a row binary tree structure
struct Row_tr {
 float e;
 int j;
 struct Row_tr *left, *right;
};

// Another binary tree holds the roots of every row
struct Header {
 Row_tr *diag;
 int i, n;                                  // row is I, row elements are n
 struct Header *top, *bottom;
};

// Les is the linear equation system class
class Les {
 int neq, iband;
 bool banded, reduced;
 struct Header *heads;
 
public:
 Les(int ne, int ib);                       // constructor functions
 ~Les();                                    // destructor functions
 Header *LocHeader(int nd1);                // locate header to row nd1
 Row_tr *Loc(int nd1, int nd);              // locate element nd2 of row nd1
 float Get(int nd1, int nd2);               // get element at nd1, nd2
 void Update(int nd1, int nd2, float newx); // change value at nd1, nd2
 void AddElem(float x, int nn1, int nn2);   // create value at nn1, nn2
 void DelRow(Row_tr *p);                    // delete row tree pointed to by p
 void DelMatrix();                          // delete system from memory
 Row_tr *RowsTree(int nel);                 // build a balanced row tree initialized to 0.0 and 0
 void Balance(int nd1);                     // balance the nd1-th row tree 
 int GetLevel(int I);                       // get the level of the I-th tree
 long Count();                              // get the total number of elements
 void Bansol(task duty);                    // solve the symmetric linear equations system

 protected:
 Header *HeadersTree(int ne);               // builds a balanced tree of ne header pointers
 void InitTree(Header *h);                  // initializes the row number on each header structure
 void newelem(float xx, int ii, Row_tr **p);// add elem to col i of row tree *p
 void DelHeads(Header (**q));               // deletes the roots of rows
 Row_tr *GetNext(int nd1);                  // get pointer to next elem of nd1-th row tree
 void CopyTree(Row_tr *p, int nd1);         // copy tree to a balanced new tree
 void level(Row_tr *a, int l);              // get the level of the a-tree
 void countelem(Row_tr *p, long *n);        // count of row elements pointed to by p
 void countrow(Header *q, long *n);         // count of system rows
};
