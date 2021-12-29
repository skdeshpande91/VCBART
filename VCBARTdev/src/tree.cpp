#include <string>
#include <vector>
#include <map>
#include "tree.h"

//--------------------------------------------------
// constructors

tree::tree(){
  mu = 0.0;
  
  rule.binary = false;
  rule.v = 0;
  rule.c = 0;

  p = 0;
  l = 0;
  r = 0;
}

// for the destructor: cut back to one node
void tree::to_null()
{
  size_t ts = get_treesize();
  while(ts > 1){
    npv nv;
    get_nogs(nv); // get nodes with no grandchildren
    for(size_t i = 0; i < nv.size(); i++){
      delete nv[i]->l;
      delete nv[i]->r;
      nv[i]->l=0;
      nv[i]->r=0;
    }
    ts = get_treesize();
  }
  mu = 0.0;
  rule.clear();
  
  p = 0;
  l = 0;
  r = 0;
}

// print
void tree::print(bool pc) const // pc is flag to print children
{
  size_t id = get_nid();
  size_t pid; // id of the parent node
  if(!p) pid = 0; // parent of this is the top node, which we'll assign a pid of 0
  else pid = p->get_nid();
  
  if(pc && (get_ntype() == 't')) Rcpp::Rcout << "tree size:" << get_treesize() << std::endl;
  Rcpp::Rcout << "id." << id;
  if(get_ntype() == 'b') Rcpp::Rcout << "  mu: " << mu << std::endl;
  else { // internal node or nnog or top node
    if(rule.binary){
      // splitting on a binary variable
      Rcpp::Rcout << "  split on binary variable Z" << rule.v+1 << std::endl;
    } else if(!rule.binary){
      Rcpp::Rcout << "  split on continuous variable Z" << rule.v+1 << " at" << c << std::endl;
    } else if(!p){
      // this is the top node
      Rcpp::Rcout << " tree contains only the root node" << std::endl;
    }
  }
  
  if(pc){
    if(l){
      l->print(pc);
      r->print(pc);
    }
  }
}



//--------------------------------------------------
//operators

// this overloads = if we say tree1 = tree2

tree& tree::operator=(const tree& rhs)
{
   if(&rhs != this) {
      to_null(); //kill left hand side (this)
      cp(this,&rhs); //copy right hand side to left hand side
   }
   return *this;
}

// tree-level gets

// get a vector of pointers to the bottom nodes
void tree::get_bots(npv& bv)
{
  if(l) { //have children
    l->get_bots(bv);
    r->get_bots(bv);
  } else bv.push_back(this);
}

//get a vector of pointers to the no grandchildren nodes
void tree::get_nogs(npv& nv)
{
  if(l) { //have children
    if((l->l) || (r->l)) {  //have grandchildren
      if(l->l) l->get_nogs(nv);
      if(r->l) r->get_nogs(nv);
    } else nv.push_back(this);
  }
}

bool tree::is_nog() const
{
  bool isnog=true;
  if(l) {
    if(l->l || r->l) isnog=false; //one of the children has children.
  } else isnog=false; //no children
  return isnog;
}

//get a vector of pointers to *ALL* nodes
void tree::get_nodes(npv& v)
{
  v.push_back(this);
  if(l) {
    l->get_nodes(v);
    r->get_nodes(v);
  }
}

void tree::get_nodes(cnpv& v)  const
{
  v.push_back(this);
  if(l) {
    l->get_nodes(v);
    r->get_nodes(v);
  }
}


// get the size of the tree (i.e. number of nodes)
int tree::get_treesize() const
{
   if(!l) return 1;  //if bottom node, tree size is 1
   else return (1+l->get_treesize()+r->get_treesize());
}
//number of nodes with no
int tree::get_nnogs() const
{
  if(!l) return 0; // this is a bottom node
  if(l->l || r->l) return(l->get_nnogs() + r->get_nnogs()); // this has at least one grandchild
  else return 1; // this is a nog node
}
// count number of bottom nodes
int tree::get_nbots() const
{
  if(!l) return 1; // this is a bottom node
  else return(l->get_nbots() + r->get_nbots());
}

// get depth of tree
int tree::get_depth() const
{
   if(!p) return 0; //no parents
   else return (1+p->get_depth());
}

// get the node id
int tree::get_nid() const
//recursion up the tree
{
   if(!p) return 1; //if you don't have a parent, you are the top
   if(this==p->l) return 2*(p->get_nid()); //if you are a left child
   else return 2*(p->get_nid())+1; //else you are a right child
}

// get the node type
char tree::get_ntype() const
{
   //t:top, b:bottom, n:no grandchildren, i:internal
   if(!p) return 't';
   if(!l) return 'b';
   if(!(l->l) && !(r->l)) return 'n';
   return 'i';
}

tree::tree_p tree::get_ptr(int nid)
{
   if(this->get_nid() == nid) return this; //found it
   if(!l) return 0; //no children, did not find it
   tree_p lp = l->get_ptr(nid);
   if(lp) return lp; //found on left
   tree_p rp = r->get_ptr(nid);
   if(rp) return rp; //found on right
   return 0; //never found it
}

tree::tree_cp tree:get_bn(double* z)
{
  // z points to the first modifier of a particular observation
  if(!l) return this; // node has no left child, it must be a bottom node
  else{
    if(z[v] < rule.c) return l->get_bn(z);
    else if(z[v] >= rule.c) return r->get_bn(z);
    else return 0; // it's a problem if we ever return 0...
  }
}

/*
 // almost certainly don't need this
double tree::evaluate(double* z){
  tree::tree_cp bn = get_bn(z);
  if(bn == 0) return std::nan(""); // no valid bottom node, return nan
  else return bn->get_mu();
}
*/



// birth
void tree::birth(int nid, rule_t rule)
{
  tree_p np = get_ptr(nid); // get the pointer to the node being split
  if(!np){
    Rcpp::Rcout << "Trying birth at node w/ id " << nid << " but pointer is invalid!" << std::endl;
    Rcpp::stop("[birth]: invalid pointer to node!");
  }
  if(np->l) Rcpp::stop("[birth]: trying to split a node that already has children!");
  
  tree_p l = new tree; // initialize the new tree that will be the left child of nid
  l->mu = 0.0; // we will overwrite this later
  l->rule.clear();

  tree_p r = new tree; // initialize the new tree that will be the left child of nid
  r->mu = 0.0; // we will overwrite this later
  r->rule.clear();
  
  np->l = l;
  np->r = r;
  
  np->rule.binary = rule.binary;
  np->rule.v = rule.v;
  np->rule.c = rule.c;
  
  np->mu = 0.0; // we will overwrite this value
  
  l->p = np;
  r->p = np;
  
}

// perform a death
void tree::death(int nid)
{
  tree_p nb = get_ptr(nid);
  if(!nb) Rcpp::stop("[death]: missing pointer for nid!");
  
  if(nb->is_nog()){
    delete nb->l;
    delete nb->r;
    
    nb->l = 0; // nb now has no children so set corresponding pointers to 0
    nb->r = 0; // nb now has no children so set corresponding pointers to 0
    nb->mu = 0.0; // this will be over-written when we update the mu's
    // reset the rule values
    nb->rule.clear();
  } else Rcpp::stop("[death]: cannot perform death move on a node with grandchildren.");
}

void tree::get_rg_cont(int &v, double& c_lower, double &c_upper){
  // starting at current node, recurse up the tree checking if parents split on v-th continuous predictor
  if(p){
    if(!p->rule.binary && p->rule.v == v){
      if(this == p->l){
        // my parent splits on v, so any observation in this should be less than p->rule.c
        // this establishes an upper bound on valid cutpoints for v in this
        // as we recurse up the tree, we don't want to overwrite a stricter region
        // if p->rule.c is somehow smaller than p->rule.c then we will overwrite it
        if(p->rule.c < c_upper) c_upper = p->rule.c;
      } else(if this == p->r){
        // my parent splits on v, so any observations in this should be greater than or equal to p->rule.c
        // this establishes a lower bound on valid cutpoints for v in this:
        if(p->rule.c > c_lower) p->rule.c;
      } else Rcpp::stop("[get_rg_cont]: this was not equal to either left or right child of parent!");
    } // closes loop checking that (i) p splits on continuous variable v
    p->get_rg_cont(v, c_lower, c_upper);
  }
}

/*
 not needed at the moment
void tree::get_rg_cat(int &v, std::set<int> &levels){
// recruse up the tree:
//    1. we will initialize levels to be the set of all levels for x_cat[v]
//    2. if this has a parent that splits on v, check if this is p->l or p->r.
//        * if this is the left child of a node that splits on v, then replace levels with p->rule.l_vals & set recurse = false
//
  bool recurse = true;
  if(p){
    // this has a parent.
    if(p->rule.is_cat && p->rule.v_cat == v){
      // parent of this splits on v
      if(this == p->l){
        levels.clear();
        for(set_it it = p->rule.l_vals.begin(); it != p->rule.l_vals.end(); ++it) levels.insert(*it);
        recurse = false; // no need to continue recursing up the tree
      } else if(this == p->r){
        levels.clear();
        for(set_it it = p->rule.r_vals.begin(); it != p->rule.r_vals.end(); ++it) levels.insert(*it);
        recurse = false; // no need to continue recursing up the tree
      } else Rcpp::stop("[get_rg_cat]: this was not equal to either left or right child of its parent!");
    }
    if(recurse) p->get_rg_cat(v, levels);
  }
}
*/
//private functions

//copy tree o to tree n
void tree::cp(tree_p n, tree_cp o)
//assume n has no children (so we don't have to kill them)
//recursion down
{
  if(n->l) Rcpp::stop("[cp]:tree n has children.");
  // if we haven't stopped by now, it's valid to continue to copying

  n->mu = o->mu;
  // not 100% sure if it's valid to do n->rule = o->rule
  // but better to be safe and deliberately copy all members
  n->rule.is_cat = o->rule.is_cat;
  n->rule.proj_vec = o->rule.proj_vec;
  n->rule.c = o->rule.c;
  n->rule.v_cat = o->rule.v_cat;
  n->rule.l_vals = o->rule.l_vals;
  n->rule.r_vals = o->rule.r_vals;
  
  if(o->l){
    // if o has children
    n->l = new tree; // create new tree for n's left child
    (n->l)->p = n; // assign parent of n's left child as n
    cp(n->l, o->l); // recurse for left child
    n->r = new tree;
    (n->r)->p = n;
    cp(n->r, o->r);
  }
}
/* [skd]: I don't think we ever use birthp and deathp*/
/*
void tree::birthp(tree_p np,size_t v, size_t c, double ml, double mr)
{
   tree_p l = new tree;
   l->mu=ml;
   tree_p r = new tree;
   r->mu=mr;
   np->l=l;
   np->r=r;
   np->v = v; np->c=c;
   l->p = np;
   r->p = np;
}
//--------------------
//kill children of  nog node *nb
void tree::deathp(tree_p nb, double mu)
{
   delete nb->l;
   delete nb->r;
   nb->l=0;
   nb->r=0;
   nb->v=0;
   nb->c=0;
   nb->mu=mu;
}
//--------------------
*/


/*
 input/output operators
 std::ostream& operator<<(std::ostream& os, const tree& t)
 {
    tree::cnpv nds;
    t.getnodes(nds);
    os << nds.size() << endl;
    for(size_t i=0;i<nds.size();i++) {
       os << nds[i]->nid() << " ";
       os << nds[i]->getv() << " ";
       os << nds[i]->getc() << " ";
       os << nds[i]->getm() << endl;
    }
    return os;
 }
 
 //input operator
 std::istream& operator>>(std::istream& is, tree& t)
 {
    size_t tid,pid; //tid: id of current node, pid: parent's id
    std::map<size_t,tree::tree_p> pts;  //pointers to nodes indexed by node id
    size_t nn; //number of nodes

    t.tonull(); // obliterate old tree (if there)

    //read number of nodes----------
    is >> nn;
    if(!is) {
       //cout << ">> error: unable to read number of nodes" << endl;
       return is;
    }

    //read in vector of node information----------
    std::vector<node_info> nv(nn);
    for(size_t i=0;i!=nn;i++) {
       is >> nv[i].id >> nv[i].v >> nv[i].c >> nv[i].m;
       if(!is) {
          //cout << ">> error: unable to read node info, on node  " << i+1 << endl;
          return is;
       }
    }

    //first node has to be the top one
    pts[1] = &t; //careful! this is not the first pts, it is pointer of id 1.
    t.setv(nv[0].v); t.setc(nv[0].c); t.setm(nv[0].m);
    t.p=0;

    //now loop through the rest of the nodes knowing parent is already there.
    for(size_t i=1;i!=nv.size();i++) {
       tree::tree_p np = new tree;
       np->v = nv[i].v; np->c=nv[i].c; np->mu=nv[i].m;
       tid = nv[i].id;
       pts[tid] = np;
       pid = tid/2;
       // set pointers
       if(tid % 2 == 0) { //left child has even id
          pts[pid]->l = np;
       } else {
          pts[pid]->r = np;
       }
       np->p = pts[pid];
    }
    return is;
 }
 std::ostream& operator<<(std::ostream& os, const xinfo& xi)
 {
   os << xi.size() << endl;
   for(size_t i=0;i<xi.size();i++) {
     os << xi[i].size() << endl;
     for(size_t j=0;j<xi[i].size();j++)
       os << xi[i][j]<<endl;
     os << endl;
   }

   return os;
 }
 std::istream& operator>>(std::istream& is, xinfo& xi)
 {
   size_t xin;
   size_t vecdn;

   xi.resize(0); // reset old xinfo (if there)

   is >> xin;
   if(!is) {
     //cout << ">> error: unable to read size of xinfo" << endl;
     return is;
   }

   std::vector<double> vec_d;
   double vecdelem;

   for(size_t i=0;i<xin;i++) {
     is >> vecdn;
     for(size_t j=0;j<vecdn;j++) {
       is >> vecdelem;
       vec_d.push_back(vecdelem);
     }
     xi.push_back(vec_d);
     vec_d.resize(0);
   }

   return is;
 }

 */
