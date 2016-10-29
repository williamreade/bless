 // High efficiency linear equations system solver module
 // Author: William Read, Universidad Nacional Pedro Henríquez Ureña, Santo Domingo.
 // May be applied, copied or reproduced, giving credit to the author.
 //
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include "bless.h"

using namespace std;

Header *CurrHead;                          // Current head
Row_tr *q;                                 // q and
int ni;                                    // ni are global for bless

// Builds a balanced tree of headers; returns pointer to the root 

Header *Les::HeadersTree(int ne)
{
  if (ne == 0) return NULL;
  else
  { int ni = ne / 2;
  int nd = ne - ni - 1;
  Header *pij = new Header;
  pij->diag = NULL;
  pij->i = 0;
  pij->n = 0;
  pij->top = HeadersTree(ni);
  pij->bottom = HeadersTree(nd);
  return pij;
  }
}

// initializes row numbers of headers structures

void Les::InitTree(Header *h)
{
  if ( h != NULL) {
  InitTree(h->top);
  h->i = ni; ++ni;
  InitTree(h->bottom);
  }
}

// the constructor function

Les::Les(int ne, int ib)
{
  neq = ne;
  iband = ib;
  banded = (neq > iband);
  heads = HeadersTree(neq);
  ni = 1; InitTree(heads);
}

// the destructor function

Les::~Les()
{
  DelMatrix();
}

// the locate header function

Header *Les::LocHeader(int nd1)
{
  Header *p = heads;
  bool found = false;
  while ((p != NULL) && (found == false)) {
    if (p->i == nd1) found = true;
    else
    if (p->i > nd1) p = p->top;
    else
    if (p->i < nd1) p = p->bottom;
  };
return p;
}

// the locate element function

Row_tr *Les::Loc(int nd1, int nd2)
{
  Row_tr *p = LocHeader(nd1)->diag;
  bool found = false;
  while ((p != NULL) && (found == false)) {
    if (p->j == nd2) found = true;
    else
    if (p->j > nd2) p = p->left;
    else
    if (p->j < nd2) p = p->right;};
  return p;
} 

// get element function

float Les::Get(int nd1, int nd2)
{
  Row_tr *p = Loc(nd1, nd2);
  if (p != NULL) return p->e;
  else
  cout << "\nThe element "<< nd1 <<", "<< nd2 <<" does not exist\n";
}

// update element

void Les::Update(int nd1, int nd2, float newxx)
{
  Row_tr *p = Loc(nd1, nd2);
  if (p != NULL) {
  p->e = newxx; return;
  }
  else
  cout << "\nThe element "<< nd1 <<", "<< nd2 <<" does not exist\n";
}

// add element at nn2 to row pointd to by p

void Les::newelem( float xx, int jj, Row_tr **p)
{
  if ((*p) == NULL) {
  (*p) = new Row_tr;
  (*p)->e = xx;
  (*p)->j = jj;
  (*p)->left = NULL;
  (*p)->right = NULL;
  }
  else
  if ((*p)->j > jj) newelem(xx, jj, &(*p)->left);
  else
  if ((*p)->j < jj) newelem(xx, jj, &(*p)->right);
}

// add element to matrix at nn1 nn2

void Les::AddElem(float x, int nn1, int nn2)
{
  CurrHead = LocHeader(nn1);
  newelem(x, nn2, &CurrHead->diag);
  ++CurrHead->n;
}

// delete the row tree pointed to by p

void Les::DelRow(Row_tr *q)
{
  if (q != NULL) {
  DelRow(q->left);DelRow(q->right);
  if ((q->left == NULL) && (q->right == NULL)) {
    delete q;
    q = NULL;
    }
  }
} 

// deletes the roots of rows

void Les::DelHeads(Header (**q))
{
  if ((*q) != NULL) {
  DelHeads(&(*q)->top);
  DelRow((*q)->diag);
  DelHeads(&(*q)->bottom);
  if (((*q)->top == NULL) && ((*q)->bottom == NULL)) {
    delete (*q);
    (*q) = NULL;
    }
  }
}

// deletes system from memory

void Les::DelMatrix() {
  DelHeads(&heads);
}

// builds a balanced row tree initialized to 0.0 and 0

Row_tr *Les::RowsTree(int nel) {
  if (nel != 0) {
  int ni = nel / 2;
  int nd = nel -ni -1;
  Row_tr *pij = new Row_tr;
  pij->e = 0.0;
  pij->j =0;
  pij->left = RowsTree(ni);
  pij->right = RowsTree(nd);
  return pij;
  }
  else
  return NULL;
}

// get the pointer to next elem of nd1-th row tree

Row_tr *Les::GetNext(int nd1) {
  Row_tr *ptst = NULL;
  while((ni < (iband+1)) && (ptst == NULL)) {
  ++ni;
  ptst = Loc(nd1, ni);
  if (ptst->j == ni) return ptst;
  }
}

// copy tree to a balanced new one

void Les::CopyTree(Row_tr *p, int nd1) {
  if (p != NULL) {
  CopyTree(p->left, nd1);
  q = GetNext(nd1);
  p->e = q->e;
  p->j = q->j;
  CopyTree(p->right, nd1);
  }
}

// balance the nd1-th row tree

void Les::Balance(int nd1) {
  CurrHead = LocHeader(nd1);
  Row_tr *qbal = RowsTree(CurrHead->n);
  ni = 0;
  CopyTree(qbal, nd1);
  DelRow(CurrHead->diag);
  CurrHead->diag = qbal;
}

// get the level of the a-tree

void Les::level(Row_tr *a, int l) {
  if (a != NULL) {
  level(a->left, l++);
  if (l > ni) ni = l;
  level(a->right, l++);
  }
}

// get the level of the i-th tree

int Les::GetLevel(int i) {
  Header *p = LocHeader(i);
  ni = 0;
  level(p->diag, ni);
  return ni;
}

// count row elements

void Les::countelem(Row_tr *p, long *n) {
  if (p != NULL) {
  countelem(p->left, &(*n));
  (*n)++;
  countelem(p->right, &(*n));
  }
}

// count of systems row

void Les::countrow(Header *q, long *n) {
  if (q != NULL) {
  countrow(q->top, &(*n));
  countelem(q->diag, &(*n));
  countrow(q->bottom, &(*n));}
}

// count elements of system

long Les::Count() {
  long number = 0;
  countrow(heads, &number);
  return number;
}

// minor of two integers

int minor(int ii, int kk) {
  if (ii < kk) return ii;
  else return kk;
}

// Banded symmetric equation system solver

void Les::Bansol(task duty) {
  int m, mr, i, j;
  float pivot, cp;
  Row_tr *a, *a1, *r, *r1;
  
// reduce system to triangle matrix

  int nrs = neq - 1;
  if (duty == trian) {
    for (int nn = 1; nn <= nrs; nn++) {
      m = nn -1;
      mr = minor(iband, neq - m);
      pivot = Get(nn, 1);
      for (int l = 2; l <= mr; l++) {
        a = Loc(nn, l);
        if (a == NULL) cp = 0.0;
        else cp = a->e / pivot;
        i = m + l;
        j = 0;
        for (int k = l; k <= mr; k++) {
          j+= 1;
          a1 = Loc(nn, k);
          if (a1 != NULL) {
              a = Loc(i, j);
              if (a == NULL) {
                AddElem(0.0, i, j);
                a = Loc(i, j);
              }
              a->e -= cp * a1->e;
            }
          }; 
          a = Loc(nn, l);
          a->e = cp;
        }
     }; 
     reduced = true;
     return;
   }

// solve reduced system for the right hand side – forward reduction

   if (duty == reslv) {
    for (int nn = 1; nn <= nrs; nn++) {
      m = nn-1;
      mr = minor(iband, neq - m);
      r = Loc(nn, iband + 1);
      a = Loc(nn, 1);
      cp = r->e;
      r->e = cp / a->e;
      for (int l = 2; l <= mr; l++) {
        i = m + l;
        r = Loc(i, iband + 1);
        a = Loc(nn, l);
        if (a != NULL) r->e -= a->e * cp;
      }
    };
    r = Loc(neq, iband +1);
    a = Loc(neq, 1);
    r->e /= a->e;
  
// back substitution

   for (i = 1; i <= nrs; i++) {
      int nn = neq - i;
      m = nn -1;
      mr = minor(iband, neq - m);
      for (int k = 2; k <= mr; k++) {
       j = m + k;
       r = Loc(nn, iband +1);
       a = Loc(nn, k);
        if (q != NULL) {
         r1 = Loc(j, iband + 1);
         r->e -= a->e * r1->e;
       }
     }
    }
  }
}
