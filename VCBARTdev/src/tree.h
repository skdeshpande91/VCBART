#ifndef GUARD_tree_h
#define GUARD_tree_h


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


// new structure for the decision rule



class tree {
public:
   //------------------------------
   //friends
  //friend std::istream& operator>>(std::istream&, tree&);
  
  //friend void grow_tree(tree &t, int &accept, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen, bool debug);
  //friend void prune_tree(tree &t, int &accept, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen, bool debug);
  //friend void update_tree(tree &t, int &accept, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen, bool debug);
  
  //------------------------------
  //typedefs
  typedef tree* tree_p;
  typedef const tree* tree_cp;
  typedef std::vector<tree_p> npv; //Node Pointer Vector
  typedef npv::iterator npv_it; // iterator for node pointer vector
  typedef std::vector<tree_cp> cnpv; //const Node Pointer Vector
  typedef cnpv::iterator cnpv_it; // iterator for const node pointer vector

  //------------------------------
  //tree constructors, destructors
  tree();
  ~tree() {to_null();}

  //------------------------------
  //operators
  tree& operator=(const tree&);

  // tree-level functions
  void to_null(); //like a "clear", null tree has just one node
  void print(bool pc = true) const;

  //node-level sets
  void set_mu(double mu) {this->mu=mu;}
  
  //node-level gets
  double get_mu() const {return mu;}
  
  bool get_binary() const {return rule.binary;} // find out if we're doing a split on a binary variable
  int get_v() const{return rule.v;}
  double get_cutpoint() const{return rule.c;}
  
  tree_p get_p() const {return p;}  //should this be tree_cp?
  tree_p get_l() const {return l;}
  tree_p get_r() const {return r;}
  
  // tree-level gets
  void get_bots(npv& bv); // get bottom nodes
  void get_nogs(npv& nv); //get nog nodes (no granchildren)
  bool is_nog() const; // is 

  void get_nodes(npv& v); //get vector of all nodes
  void get_nodes(cnpv& v) const; //get all nodes
  
  int get_treesize() const; // get number of nodes in tree
  int get_nnogs() const; // get number of nog nodes (nodes with no grandchildren)
  int get_nbots() const; // get number of bottom nodes
  
  int get_depth() const; //depth of a node
  int get_nid() const;   //node id
  char get_ntype() const;   //t:top;b:bottom;n:nog;i:interior, carefull a t can be bot
  tree_p get_ptr(int nid); //get node pointer from node id, 0 if not there.
  //tree_cp get_ptr(int nid); // node pointer
  tree_cp get_bn(double* z); // get the pointer corresponding to the bottom node containing x
  
  //double evaluate(double* z); // evaluate tree at a single point
  

  //birth death using nid----------
  void birth(int nid, rule_t rule);
  void death(int nid);

  void get_rg_cont(int &v, double& c_lower, double &c_upper)
  //void get_rg_cat(int &v, std::set<int> &levels); // what are the available levels of a categorical variable at some node
  
private:
   //------------------------------
   //parameter for node
  double mu;
  rule_t rule;

  //------------------------------
  //tree structure
  tree_p p; //parent
  tree_p l; //left child
  tree_p r; //right child
  //------------------------------
  //utiity functions
  
  // [skd]: unlikely that we actually use these utility functions
  
  void cp(tree_p n,  tree_cp o); //copy tree o to n
  //void birthp(tree_p np, rule_t rule);
  //void deathp(tree_p nb, double mu); //kill children of nog node nb
};
//std::istream& operator>>(std::istream&, tree&);
//std::ostream& operator<<(std::ostream&, const tree&);
//std::istream& operator>>(std::istream&, xinfo&);
//std::ostream& operator<<(std::ostream&, const xinfo&);


#endif
