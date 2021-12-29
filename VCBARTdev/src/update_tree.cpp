//
//  update_tree.cpp
//  
//
//  Created by Sameer Deshpande on 12/29/21.
//

#include "update_tree.hpp"


void grow_tree( const suff_stat &orig_suff_stat, suff_stat &new_suff_stat,){
  
  
  tree::npv bnv;
  t.get_bots(bnv); // get all of the bottom nodes
  
  // pick one of the bottom nodes
  tree::tree_cp nx = ;
  
  // once we have picked a node
  int nx_id = nx->get_nid(); // get the id of the node we are going to try to split
  
  int nxl_id = 2*nx_id; // id of the proposed left child
  int nxr_id = 2*nx_id + 1; // id of the proposed right child
  
  // for a grow move
  
  
  
  
  
  // once we draw a new value of mu
  // we need to do something like
  tree::tree_p nxl = t.get_ptr(nxl_id);
  tree::tree_p nxr = t.get_ptr(nxr_id);
  
  nxl->setmu();
  nxr->setmu();
  
  
}
