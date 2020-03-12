#ifndef GUARD_tree_h
#define GUARD_tree_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <vector>

#include "info.h"
#include "rng.h"

/*
Class is a bit confusing because sometimes you
want to think about an instance as a node and
sometimes as the whole tree including the node
and all its children.
Fundamentally, we have a tree not a node.
*/

/*
three ways to access a node:
(i) node id: is the integer assigned by the node numbering system 
     assuming they are all there
(ii) node ind: is the index into the array of 
   node(tree) pointers returned by getnodes or getbots or getnogs
   which means you go left to right across the bottom of the tree
(iii) by its pointer (should only "give out" const pointers)
*/

//info contained in a node, used by input operator
struct node_info {
   std::size_t id; //node id
   std::size_t v;  //variable
   std::size_t c;  //cut point
   double m;       //mu
};

//xinfo xi, then xi[v][c] is the c^{th} cutpoint for variable v.
// left if x[v] < xi[v][c]
typedef std::vector<double> vec_d; //double vector
typedef std::vector<vec_d> xinfo; //vector of vectors, will be split rules

class tree {
public:
   //------------------------------
   //friends
  friend std::istream& operator>>(std::istream&, tree&);
  friend void bd_ind(tree &x, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen); // for independent error correlation structures
  friend void bd_cs(tree &x, double &rho, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen); // for CS error correlation structure
  friend void bd_ar(tree &x, double &rho, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen); // for AR error correlation structure
  
  //friend double bd_ar(tree &x, std::vector<arma::mat> Omega_vec, std::vector<double> &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);
  //friend double bd_ar(tree &x, std::vector<arma::mat> Omega_vec, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);
  //friend double bd_ind(tree &x, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);
  //friend void bd_alt(tree &x, std::vector<arma::mat> Omega_vec, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen, size_t min_leaf);
  //friend void bd_cs(tree &x, double &rho, double &sigma, xinfo &xi, data_info &di, tree_prior_info &tree_pi, RNG &gen);
 
  
   //------------------------------
   //typedefs
   typedef tree* tree_p;
   typedef const tree* tree_cp;
   typedef std::vector<tree_p> npv; //Node Pointer Vector
   typedef std::vector<tree_cp> cnpv; //const Node Pointer Vector

   //------------------------------
   //tree constructors, destructors
   tree();
   tree(const tree&);
   tree(double);
   ~tree() {tonull();}

   //------------------------------
   //operators
   tree& operator=(const tree&);

   //------------------------------
   //node access
   //you are freely allowed to change mu, v, c
   //set----------
   void setm(double mu) {this->mu=mu;}
   void setv(size_t v) {this->v = v;}
   void setc(size_t c) {this->c = c;}
   //get----------
   double getm() const {return mu;} 
   size_t getv() const {return v;}
   size_t getc() const {return c;}
   tree_p getp() const {return p;}  //should this be tree_cp? 
   tree_p getl() const {return l;}
   tree_p getr() const {return r;}

   //------------------------------
   //tree functions
   //stuff about the tree----------
   size_t treesize() const; //number of nodes in tree
   size_t nnogs() const;    //number of nog nodes (no grandchildren nodes)
   size_t nbots() const;    //number of bottom nodes
   void pr(bool pc=true) const; //to screen, pc is "print children"
  void print(bool pc = true) const;
  
  
   //birth death using nid----------
   bool birth(size_t nid,size_t v, size_t c, double ml, double mr); 
   bool death(size_t nid, double mu); 
   //vectors of node pointers----------
   void getbots(npv& bv);         //get bottom nodes
   void getnogs(npv& nv);         //get nog nodes (no granchildren)
   void getnodes(npv& v);         //get vector of all nodes
   void getnodes(cnpv& v) const;  //get all nodes
   //find node from x and region for var----------
   tree_cp bn(double *x,xinfo& xi);  //find bottom node for x
   void rg(size_t v, int* L, int* U) const; //find region [L,U] for var v.
   size_t nuse(size_t v); //how many times var v is used in a rule.
   void varsplits(std::set<size_t> &splits, size_t v); //splitting values for var v.
   //------------------------------
   //node functions
   size_t depth() const; //depth of a node
   size_t nid() const;   //node id
   char ntype() const;   //t:top;b:bottom;n:nog;i:interior, carefull a t can be bot
   tree_p getptr(size_t nid); //get node pointer from node id, 0 if not there.
   bool isnog() const;
   void tonull(); //like a "clear", null tree has just one node

private:
   //------------------------------
   //parameter for node
   double mu; 
   //------------------------------
   //rule: left if x[v] < xinfo[v][c]
   size_t v; 
   size_t c; 
   //------------------------------
   //tree structure
   tree_p p; //parent
   tree_p l; //left child
   tree_p r; //right child
   //------------------------------
   //utiity functions
   void cp(tree_p n,  tree_cp o); //copy tree o to n
   void birthp(tree_p np,size_t v, size_t c, double ml, double mr); 
   void deathp(tree_p nb, double mu); //kill children of nog node nb 
};
std::istream& operator>>(std::istream&, tree&);
std::ostream& operator<<(std::ostream&, const tree&);
std::istream& operator>>(std::istream&, xinfo&);
std::ostream& operator<<(std::ostream&, const xinfo&);


#endif
