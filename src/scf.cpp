#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <string>
#include <unordered_map>
#include "scf.h"
using namespace std;
using namespace arma;
SCF::SCF(){
  //Read in the object name like H2O
  std::string path_data = "/home/xiuyiqin/GF3ES/Test_SCF/data/";
  ifstream iptcoor(path_data+"inputfile.com");
  //In the input file never try to turn on the ERI symmetry, it will cause errors
  std::string data;
  string judge1 = "Obj=";
  string judge2 = "OCC=";
  while(!iptcoor.eof()){
      iptcoor >> data;
        if (data.compare(judge1)==0){
            iptcoor>>data;
            objname = data;
        }
        if (data.compare(judge2)==0){
            iptcoor>>data;
            _numoccp = std::stoi(data,nullptr,10);
        }
  } 
  if (_numoccp == 0){
      printf("Please put number of occupied orbits in inputfile.com");
      exit(EXIT_FAILURE);
  }
  //Read in the 3 input files
  std::string overlapfile = "_overlap.dat";
  std::string coreHfile = "_coreH.dat";
  std::string erifile = "_eri.dat";
  std::string enucfile = "_enuc_Nb.dat";
  ifstream iptovlap(path_data+objname+overlapfile);
  ifstream iptcoreH(path_data+objname+coreHfile);
  ifstream ipteri(path_data+objname+erifile);
  ifstream iptenuc(path_data+objname+enucfile);
  
  while(!iptenuc.eof()){
  iptenuc>>_nucrepul;
  iptenuc>>_numorbit;
  }
  //to set up ioff:
  _ioff.resize(pow(_numorbit,4));
  _ioff[0] = 0;
  for(int i=1; i < pow(_numorbit,4); i++)
  _ioff[i] = _ioff[i-1] + i;
  //set up matrix
  _ovlap = new arma::mat(_numorbit,_numorbit);
  _coreHam = new arma::mat(_numorbit,_numorbit);
  _Pini = new arma::mat(_numorbit,_numorbit);
  _Pnext = new arma::mat(_numorbit,_numorbit);
  _Fini = new arma::mat(_numorbit,_numorbit);
  _Fnext = new arma::mat(_numorbit,_numorbit);
  int i = 0, j = 0,k = 0, l = 0;
  double temp,temp2;
  printf("T3");
  /*
  to input overlap integral
  and core Hamiltonian
  and nuclear repulsive energy and No of orbits 
  */
   // overlap and core Hamiltonian is same size so I put them in same while
  while(!iptovlap.eof()){
  iptovlap >> i >> j >> temp;
  iptcoreH >> i >> j >> temp2;
  (*_ovlap)(i -1, j -1) = temp;
  (*_ovlap)(j -1, i -1 ) = temp;
  (*_coreHam)(i -1, j -1) = temp2;
  (*_coreHam)(j -1, i -1 ) = temp2;
  }
  /*two electron integral */
  double temp4;
  while (!ipteri.eof()) {
    ipteri >> i >> j >> k >> l >> temp4;
    _twoelec[getijkl(i, j, k, l)] = temp4;
  }
  iptcoor.close();iptovlap.close();iptcoreH.close();
  ipteri.close();iptenuc.close();
}

int SCF::getijkl(int i, int j, int k, int l){
  int ij = 0, kl = 0, ijkl = 0;
  ij = (i > j) ? _ioff[i] + j : _ioff[j] + i;
  kl = (k > l) ? _ioff[k] + l : _ioff[l] + k;
  ijkl = (ij > kl) ? _ioff[ij] + kl : _ioff[kl] + ij;
  return ijkl;
}

SCF::~SCF(){
  delete _ovlap;
  delete _coreHam;
}

void SCF::print(arma::mat ipt){
  for (size_t i = 0; i < _numorbit; i++) {
    printf("%s\n", " ");
    for (size_t j = 0; j < _numorbit; j++) {
    printf("%20.12f, %s", (ipt)(i, j), " ");
    }
  }
}


void SCF::calculation(){
  /*diagonlize overlap Matrix*/
  printf("T9");
  arma::vec Seigval;
  arma::mat Seigvec;
  arma::eig_sym(Seigval, Seigvec, *_ovlap);
  /*to get diagonlized eigenvalue matrix
  *arma::mat eigvalmat = eigvec.t() * (*_ovlap) * eigvec;//works good
  but the following one would be easier to use*/
  arma::mat Seigvalmat = arma::diagmat(Seigval);
  // to get Orthogonalization Matrix
  //So element-wise inverse and square-root ! except ij term in matrix!

  /*arma::mat et = 1/ sqrt(abs(eigvalmat));
  *arma::mat S_ihalf = arma::eye(_numorbit,_numorbit);
  for (size_t i = 0; i < _numorbit; i++) {
      S_ihalf(i, i) = sqrt(1.0 /eigvalmat(i,i));
  } this one works but the following one would be more neat
  */
  arma::mat lmd_sqrti = arma::sqrt(arma::inv(Seigvalmat));

  _Ssqrtinv = Seigvec * lmd_sqrti * Seigvec.t();//different from website But I believe it caused by different order of eigenvalue? No it is caused by wrong input file
  /*temporaly moving on
  * get the Fock matrix for inital guess from hailtonian
  * then transfer to orthogonal basis to get Fock prime
  * then diagonlize Fock! in the orthogonal basis!
  * then transfer back to original basis which is the atomic orbitals
  */
  *_Fini = *_coreHam;
  _Finiprime = _Ssqrtinv.t() * (*_coreHam) * _Ssqrtinv;
  //digonalize Fock Matrix
  arma::vec Feigval;
  arma::mat Feigvec;
  arma::eig_sym(Feigval, Feigvec, _Finiprime);
  arma::mat Coe_orig = _Ssqrtinv * Feigvec;
  //density matrix equals to C * C.t()
  *_Pini = Coe_orig.cols(0, _numoccp -1) *  Coe_orig.cols(0, _numoccp -1).t();
  //then get the Energy
  arma::mat Eini = (*_Pini) % ((*_coreHam) + *_Fini);
  _Eini = arma::accu(Eini);
  _Eini += _nucrepul;
  printf("%f\n", _Eini);

  /*iteration!
  */
  while(true){
  /*build new Fock ! by old density matrix */
  for (size_t i = 0; i < _numorbit; i++) {
      for (size_t j = 0; j < _numorbit; j++) {
        (*_Fnext)(i,j) = (*_coreHam)(i,j);
        for (size_t k = 0; k < _numorbit; k++) {
          for (size_t l = 0; l < _numorbit; l++) {
            (*_Fnext)(i,j) += (*_Pini)(k,l) * (2 * _twoelec[getijkl(i+1,j+1,k+1,l+1)] - _twoelec[getijkl(i+1,k+1,j+1,l+1)]);
          }
        }
      }
    }
    _Finiprime = _Ssqrtinv.t() * (*_Fnext) * _Ssqrtinv;
    /* to get density matrix and get Energy again..
    */
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, _Finiprime);
    arma::mat Coe = _Ssqrtinv * eigvec;
    //density matrix equals to C * C.t()
    *_Pnext = Coe.cols(0, _numoccp -1) *  Coe.cols(0, _numoccp -1).t();
    //print(*_Pnext);
    _Enext = arma::accu((*_Pnext) % ((*_coreHam) + *_Fnext)) + _nucrepul;
    printf("%20.12f\n", _Enext);
    if(abs(_Enext - _Eini) < 0.00000000001)
    break;
    _Eini = _Enext;
    *_Fini = *_Fnext;
    *_Pini = *_Pnext;
  }
}
