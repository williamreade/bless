// Test-program: Builds a large symmetric banded linear equation system of input ne equations an input iband.
// The test system consists of 10's in the main diagonal and 1's at both sides of the main diagonal.
// The right hand side is the sum of all terms of each row.
// Hence the result must be all one's for the solution vector.
//
Les *Coef;
void InputMat(int ne, int ib) {
  float xl;
  for (int im = 1; im <= ne; im ++) {
    int mr = minor(ib, ne - im +1);
    for (int jm = 1; jm <= mr; jm++) {
      if (jm == 1)
      xl = 10.0;
      else
      xl = 1.0;
      Coef->AddElem(xl, im, jm);
    }
    Coef->Balance(im);
  };
}
  
void PrintMat(int ne, int ib) {
  for (int i = 1; i <= ne; i++) {
    int mr = minor(ib, ne - i + 1);
    for (int j = 1; j <= mr; j++) {
      float xl = Coef->Get(i, j); out << setiosflags(ios::left)
       << setw(10)
       << setprecision(6);
      cout << xl << " ";
    }; cout << "\n";
  }; return;
}

void InputRHS(int ne, int ib) {
  int n1 = ne -ib +1;
  int i3 = ib + 8;
  for (int i = 1; i <= ne; i++) {
    if (i <= ib) i3 += 1;
    if (i > n1) i3 -= 1;
    Coef->AddElem(i3*1.0, i, ib+1);
  }
}

void PrintRHS(int ne, int ib) {
  for(int i = 1; i <= ne; i++) cout << Coef->Get(i, ib+1) <<"\n";
}

main() {
  int ne, ib;
  char ch;
  cout << "\nNumber of equations ";
  cin >> ne;
  cout << "\nHalf bandwidth ";
  cin >> ib;
  Coef = new Les(ne, ib);
  InputMat(ne, ib);
  PrintMat(ne, ib);
  cout << "In total were generated "<< Coef->Count() << " elements\n";
  cout << "The maximum level of the tree is " << Coef->GetLevel(1) <<"\n";
  cin >> ch;
  InputRHS(ne, ib);
  PrintRHS(ne, ib);
  cin >> ch;
  Coef->Bansol(trian);
  PrintMat(ne, ib);
  cin >> ch;
  Coef->Bansol(reslv);
  PrintRHS(ne, ib);
  Coef->~Les();
}
