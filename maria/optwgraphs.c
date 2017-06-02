//
//  optwgraphs.c
//
//
//  Created by Maria Cameron on 6/2/15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.141592653589793
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))

#define INFTY 1.0e+200
#define N  1000000/* the maximal admissible number of states */
#define N_ARCS 2000000 /* the maximal admissible number of arcs */
#define STOP 0 // run algorithm while the main tree contains > STOP arcs
#define WRITE_FILE 'n' // 'y'('n'): do (do not) create files with optimal W graphs
#define MAINCYCLEMAX 1000000 /* the maximal number of states in the main cycle */
#define MAX_CONNECTED_COMPONENTS 1000000 /* the maximal admissible number of connected components */

struct my_arc {
    int tail;
    int head;
    double cost;
};

struct my_treenode { // tree node of a binary tree
    struct my_treenode *parent;
    struct my_treenode *leftchild;
    struct my_treenode *rightchild;
    int ind; // the index of the treenode
};

struct mylinklist {
    int index; // the index of the edge
    struct mylinklist *next;
    struct mylinklist *prev;
    struct mylinklist *first;
    struct mylinklist *last;
    int count; // the number of edges in the link list */
};

struct simplelist {
    int index;
    struct simplelist *next;
    struct simplelist *last;
};

// MAIN FUNCTIONS
int main(void);
void readinput(void);
void init(void);
void add_arc(void); // adds an arc into the optimal w-graph
void maximal_connected_components(void); // finds sizes of connected components

// auxiliary functions
void maketreenode( struct my_treenode *node,int ind );
int addtree(struct my_treenode *mother, struct my_treenode *node, int *count ); // returns index of the next root
int deltree( struct my_treenode *tree, int *count ); // returns index of the next root
int updatetree( struct my_treenode *tree, struct my_treenode *node ); // returns index of the next root
void swap( struct my_treenode *parent, struct my_treenode *child, char ch ); // swaps parent and its child in a binary tree
void update_arcs(struct my_treenode *tree,double delta_u); /* updates arcs' weights when a cycle is encountered */
void merge_linklists(struct mylinklist *list1, struct mylinklist *list2, int ind, char ch);
void make_link(int ind, char ch);
int find_cycle(int tail, int head, int newarc);
void merge_vlists(int head,int tail);
void find_optimal_wgraph(int sink);

//---- the input variables -----
int neib[N]; // the # of outgoing arcs of vertices
struct my_arc *arc;  // the set of arcs in the graph
const char *funame = "LJ14outgoing_arcs_weights.txt";
const char *fneibname = "LJ14outgoing_arcs.txt";

//---- variables for the binary tree ----
struct my_treenode *treenode; // each arc is assosiated with a tree node

//---- variables for computing optimal W-graphs --------
int k_wgraph; // the current index of the optimal W-graph
int nstates, narcs; // the actual numbers of states and arcs
int min_arc_index[N]; // min_arc_index[i] = the index of the minimal outgoing arc from i
int next_arc_index[N]; // next_arc_index[i] = the index of next minimal outgoing arc from i
int *count_vertex; // the number of nodes in the tree B_i of remaining outgoing arcs from the vertex i
int *count_main; // the number of nodes in the main tree
int root_index; // the index of the arc just removed from the top of the main tree
double delta[N]; // numbers Delta
struct mylinklist *first[N]; // pointers to the linklists of representing arcs in the connected components of the optimal W-graphs

// ---- variables attached to cycles ----
struct mylinklist *cfirst[N]; // pointers to the linklists of arcs forming cycles
int ncycles; // the total number of cycles

//---- variable attached to vertices ----
struct simplelist *vfirst[N]; // pointers to linklists of connected components
// the first vertex in the list is always the sink of the connected component


//---- variables for finding optimal W-graphs ----
struct simplelist *inlist;
struct simplelist *infirst[N]; // lists of incoming edges from the vertices that were new_arc
char used[N];

//---- output variables ----
double delta[N];
const char *foutname = "g_output.txt";
const char *fcyclesname = "cfirst.txt";
const char *farcsname = "used_arcs.txt";
int n_connected_components = 0; // the # of connected components in the optimal W-graph with the minimal number of sinks
int main_sink[MAX_CONNECTED_COMPONENTS]; // sinks of connected components
int size_con_comp[MAX_CONNECTED_COMPONENTS]; // the numbers of states in the connected components
int cc_count; // counter of states in connected components
int maincycle[N];
int nsub[N];
int nmaincycle = 0;




FILE *fg, *fc;


/***********************   r e a d     i n p u t   ************************/

void readinput() {
    int j,k,n,ned,ind,ind0,tsindex;
    FILE *fu, *fn, *fin;
    double u;
    
    printf("IN READINPUT\n");
    
    treenode = (struct my_treenode *)malloc(N_ARCS*sizeof(struct my_treenode));
    arc = (struct my_arc *)malloc(N_ARCS*sizeof(struct my_arc));
    count_vertex = (int *)malloc(N*sizeof(int));
    count_main = (int *)malloc(sizeof(int));
    
    fu = fopen(funame,"r");
    fn = fopen(fneibname,"r");
    
    if( fn == NULL || fu == NULL ) {
        printf("Cannot find input files\n");
        exit(1);
    }
    else {
        n = 0;
        ind = 0;
        while( !feof(fn) && !feof(fu) ) {
            count_vertex[n] = 0;
            fscanf(fn,"%i\t",&k);
            fscanf(fu,"%i\t",&k);
            neib[n] = k;
            if( k > 0 ) {             /* define the first outgoing arc for the vertex n */
                fscanf(fn,"%i\t",&k);
                (arc + ind) -> tail = n;
                (arc + ind) -> head = k;
                fscanf(fu,"%le\t",&u);
                (arc + ind) -> cost = u;
                maketreenode(treenode + ind,ind);
                ind0 = ind;
                ind++;
                count_vertex[n]++;
                if( neib[n] > 1 ) { /* define other outgoing arcs */
                    for( j = 1; j < neib[n]; j++ ) {
                        fscanf(fn,"%i\t",&k);
                        (arc + ind) -> tail = n;
                        (arc + ind) -> head = k;
                        fscanf(fu,"%le\t",&u);
                        (arc + ind) -> cost = u;
                        maketreenode(treenode + ind,ind);
                        ind0 = addtree(treenode + ind0,treenode + ind,count_vertex + n);
                        ind++;
                    }
                }
                min_arc_index[n] = (treenode + ind0) -> ind;
            }
            else { /* the vertex is the main sink of the corresponding connected component */
                printf("No outgoing arcs from state %i\n",n);
                main_sink[n_connected_components] = n;
                n_connected_components++;
                min_arc_index[n] = -1;
            }
            n++;
        } /* end while */
        nstates = n;
        narcs = ind;
        printf("nstates = %i\n narcs = %i\n",nstates,narcs);
        if( n_connected_components > 0 ) {
            printf("%i states have no outgoing arcs\n",n_connected_components);
        }
        fclose(fn);
        fclose(fu);
        
    }
    
}

/*************************** i n i t i a l i z a t i o n ********************************/

void init() {
    int i,ind,n0;
    
    k_wgraph = nstates;
    ncycles = 0;
    
    printf("IN INIT\n");
    // find the first vertex having outgoing arcs
    n0 = 0;
    while( neib[n0] == 0 ) n0++;
    // delete the min_arc from the tree of outgoing arcs from this vertex and add it to the main tree B
    // next_arc_index[state] = index of the outgoing arc on the top of the tree B_i
    next_arc_index[n0] = deltree(treenode + min_arc_index[n0],count_vertex + n0);
    root_index = min_arc_index[n0];
    maketreenode(treenode + min_arc_index[n0], min_arc_index[n0]);
    *count_main = 1;
    // delete min_arcs from trees of  outgoing arcs of the rest of the vertices and add them to the main tree
    for( i = n0 + 1; i < nstates; i++ ) {
        while( neib[i] == 0 ) i++; // if state has no outgoing arcs, skip it
        next_arc_index[i] = deltree(treenode + min_arc_index[i],count_vertex + i);
        maketreenode(treenode + min_arc_index[i],min_arc_index[i]);
        root_index = addtree(treenode + root_index,treenode + min_arc_index[i],count_main);
    }
    if( *count_main != nstates - n_connected_components ) {
        printf("Error in init: %i arcs in the main tree IS NOT EQUAL TO\n %i states - %i states with no outgoing arcs\n",
               *count_main,nstates,n_connected_components);
        exit(1);
    }
    
    // set up linklists of connected components of the hierarchy of optimal W-graphs
    for( i = 0; i < nstates; i++ ) {
        vfirst[i] = (struct simplelist *)malloc(sizeof(struct simplelist));
        vfirst[i] -> index = i;
        vfirst[i] -> next = NULL;
        vfirst[i] -> last = vfirst[i];
    }
}

/*******************   a d d     a r c   *********************/
void add_arc() {
    int tail, head, k, root, next_arc, new_arc, isink;
    int ncycles = 0, nnew_arcs = 0, last_sink;
    double delta_u; // the amount by which Uij should be increased
    char ch;
    FILE *fout, *fnew, *fcf, *fss;
    int stop = STOP; /* stop if the main tree contains stop arcs */
    struct mylinklist *linklist;
    struct simplelist *vlist;

    
    printf("IN ADD_ARC\n");
    fout = fopen(foutname,"w");
    fcf = fopen(fcyclesname,"w");
    fnew = fopen(farcsname,"w");
    while( (*count_main) > stop ) {
        // remove the top arc from the main tree and denote it by new_arc
        new_arc = root_index;
        root_index = deltree(treenode +root_index,count_main);
        tail = (arc + new_arc) -> tail;
        head = (arc + new_arc) -> head;
//        printf("%i\t%i\t%.10e\n",tail,head,(arc + new_arc) -> cost);
        fprintf(fnew,"%i\t%i\t%.10e\n",tail,head,(arc + new_arc) -> cost);
        nnew_arcs++;
        if( first[head]  == NULL  && first[tail] == NULL ) { // start new linklist
            ch = 0;
            make_link(new_arc,'l');
        }
        else if( first[head] != NULL && first[tail] == NULL ) { // add the arc to head's linklist
            ch = 1;
            merge_linklists(first[head],first[tail],new_arc,'l');
        }
        else if( first[tail] != NULL && first[head] == NULL ) { // add the arc to tail's linklist
            ch = 2;
            merge_linklists(first[tail],first[head],new_arc,'l');
        }
        else if( first[tail] != first[head] ) { // merge two linklists */
            ch = 3;
            merge_linklists(first[head],first[tail],new_arc,'l');
         }

        // A CYCLE IS CREATED
        
        else if( first[head] == first[tail] ) { // a cycle is created
            ch = 4;
            ncycles++;
            if( cfirst[tail]!= NULL && cfirst[head] != NULL && cfirst[tail] == cfirst[head] ) {
                printf("Error: cfirst[tail] == cfirst[head], %i --> %i\n", (arc + root_index) -> tail, (arc + root_index) ->tail);
                exit(1);
            }
            else {
                root = find_cycle(tail,head,new_arc);
            }
            // delete the top arc if its tail and head are both in the same cycle
            while( cfirst[(arc + root) -> tail] == cfirst[(arc + root) -> head] && root >= 0 ) {
                root = deltree(treenode + root,count_vertex + tail);
            }
            
            if( root >= 0 ) {
                next_arc = deltree(treenode + root,count_vertex + tail);
                linklist = cfirst[tail];
                do {
                    min_arc_index[(arc + (linklist -> index)) -> tail] = root;
                    min_arc_index[(arc + (linklist -> index)) -> head] = root;
                    next_arc_index[(arc + (linklist -> index)) -> tail] = next_arc;
                    next_arc_index[(arc + (linklist -> index)) -> head] = next_arc;
                    count_vertex[(arc + (linklist -> index)) -> tail] = count_vertex[tail];
                    count_vertex[(arc + (linklist -> index)) -> head] = count_vertex[tail];
                    linklist = linklist -> next;
                }
                while( linklist != NULL );
                maketreenode(treenode + root,root);
                // attach the new treenode to the main tree
                if( root_index >= 0 ) root_index = addtree(treenode + root_index,treenode + root,count_main);
                else root_index = addtree(treenode + root_index,NULL,count_main);
                
            }
            else {
                main_sink[n_connected_components] = tail;
                n_connected_components++;
            }
            
        } // end else if cycle is created
        
        // IF NO CYCLE IS CREATED: record data and merge connected components
        if( ch < 4 ) {
            delta[k_wgraph] = (arc + new_arc) -> cost;
            // find the number of nodes in the connected component being absorbed
            vlist = vfirst[tail] -> next;
            k = 1;
            while( vlist != NULL ) {
                vlist = vlist -> next;
                k++;
            }
            fprintf(fout,"%i\t%.10f\t%i\t%i\t%i\t",k_wgraph,delta[k_wgraph],vfirst[tail] -> index,tail,head);
            fprintf(fout,"%i\t",k);
            // merge connected components
            merge_vlists(head,tail);
            fprintf(fout,"%i\n", vfirst[tail] -> index);
            last_sink = vfirst[tail] -> index;
            k_wgraph--;
        }

        // IN ANY CASE: update list of incoming arcs
        inlist = (struct simplelist *)malloc(sizeof(struct simplelist));
        inlist -> index = new_arc;
        inlist -> next = NULL;
        inlist -> last = inlist;
        if( infirst[head] == NULL ) {
            infirst[head] = inlist;
        }
        else {
            (infirst[head] -> last) -> next = inlist;
            infirst[head] -> last = inlist;
        }
        
        // IF ALL OPTIMAL W-GRAPHS ARE FOUND: record all cycles
        if( (*count_main) == 1 ) {
            for( k = 0; k < N; k++ ) {
                if( cfirst[k] == NULL ) {
                    fprintf(fcf,"%i\t%i\n",k,k);
                }
                else {
                    fprintf(fcf,"%i\t%i\n",k,(arc + (cfirst[k] -> index)) -> head );
                }
            }
        }
        
        if( *count_main == 1 ) {
            main_sink[n_connected_components] = last_sink;
            n_connected_components++;
            maximal_connected_components();
        }
        
    } // end while
    
    fprintf(fout,"%i\t%10f\t%i\t0\t0\t%i\t1\n",k_wgraph,0.0,vfirst[head] -> index,nstates);
    fclose(fout);
    fclose(fcf);
    fclose(fnew);
    printf("the # of cycles: %i; the # of arcs used: %i\n",ncycles,nnew_arcs);
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/**********************   m a k e    t r e e     n o d e   ******************************/
void maketreenode( struct my_treenode *node,int ind ) {

    node -> parent = NULL;
    node -> leftchild = NULL;
    node -> rightchild = NULL;
    node -> ind = ind;

}


/**********************   a d d    a    n o d e    t o    t h e    t r e e   ************************/
int addtree(struct my_treenode *mother, struct my_treenode *node, int *count) {
    // addtree(...) adds binary tree hanging down from the node (the node and its children, grandchildren, etc)
    // to the binary tree containing node *mother (not necessarily its root).
    // addtree(...) returns the index of the root of the resulting tree.
    // addtree(...) adds node and its children, its grandchildren, etc, one in a time recursively.
    
    struct my_treenode *leftchild, *rightchild, *parent;
    int c, k, ind;
    
    
    if( node == NULL && mother != NULL) {
        return mother -> ind;
    }
    if( mother == NULL && node != NULL) {
        mother = node;
        return mother -> ind;
    }
    if( mother == NULL && node == NULL ) {
        return -1;
    }

    // Continue if *node is not NULL and *mother is not NULL
    // Save node's children
    leftchild = node -> leftchild;
    rightchild = node -> rightchild;
    node -> leftchild = NULL;
    node -> rightchild = NULL;
    (*count)++;

    parent = mother;
    
    /* find the root of the tree */
    while( parent -> parent != NULL ) {
        parent = parent -> parent;
    }
    
    // find the first vacancy in the tree and the path to it from the root */
    // labeling of the tree: root: 1; its children: 2 and 3; etc.
    k = 1;
    c = 0;
    while( k < (*count) ) {
        k *= 2;
        c++;
    }
    if( (*count) < k ) c--;
    c--;
    while( c >= 1 ) {
        k = (*count) >> c;
        if( k & 1 ) {
            if( parent -> rightchild == NULL ) {
                printf("Error in addtree: parent -> ind = %i, rightchild = NULL\n",parent -> ind);
                exit(1);
            }
            parent = parent -> rightchild;
        }
        else {
            if( parent -> leftchild == NULL ) {
                printf("Error in addtree: parent -> ind = %i, leftchild = NULL\n",parent -> ind);
                exit(1);
            }
            parent = parent -> leftchild;
        }
        c--;
    }
    k = (*count) >> c;
    // place the new node in the tree
    node -> parent = parent;
        if( k & 1 ) {
        parent -> rightchild = node;
    }
    else {
        parent -> leftchild = node;
    }
    
    
    /* find a right position for the added node in the tree */
    while( node -> parent != NULL &&
          ((arc + (node -> ind)) -> cost) < ((arc + ((node -> parent) -> ind)) -> cost) ) {
        if( node == (node -> parent) -> leftchild ) {
            swap(node -> parent, node, 'l');
        }
        else if( node == (node -> parent) -> rightchild ) {
            swap(node -> parent, node, 'r');
        }
        else {
            printf("Error in addtree\n");
            exit(1);
        }
    }
    
    node = mother;
    // find the root of the tree
    while( node -> parent != NULL ) {
        node = node -> parent;
    }
    ind = node -> ind;
    
    // add the children of the node to the tree
    if( leftchild != NULL )
        ind = addtree( mother, leftchild, count );
    if( rightchild != NULL )
        ind = addtree( mother, rightchild, count );
    
    
    return ind;
    
}

/*********   d e l e t e     t h e     r o o t     o f     t h e     t r e e   **********/
int deltree(struct my_treenode *tree, int *count) {
    // The last node in the binary tree is deleted from its position and moved to the top. Then its position is adjusted.
    struct my_treenode *leftchild, *rightchild, *node;
    int k, c, ind;
    double uleft, uright, u;
    
    if( (*count) < 1 ) {
        printf("Error in deltree: count = %i\n",*count);
        exit(1);
    }
    
    if( (*count) == 1 ) {
        (*count)--;
        tree = NULL;
        return -1;
    }
    
    leftchild = tree -> leftchild;
    rightchild = tree -> rightchild;
    
    node = tree;
    // find the path to the last node of the tree from the root
    k = 1;
    c = 0;
    while( k < (*count) ) {
        k *= 2;
        c++;
    }
    if( (*count) < k ) c--;
    c--;
    while( c >= 1 ) {
        k = (*count) >> c;
        if( k & 1 ) {
            node = node -> rightchild;
        }
        else {
            node = node -> leftchild;
        }
        c--;
    }
    k = (*count) >> c;
    // Delete the last node from it position.
    if( k & 1 ) {
        node = node -> rightchild;
        (node -> parent) -> rightchild = NULL;
    }
    else {
        if( node -> leftchild == NULL ) {
            printf("Error in deltree: last left:\n");
            exit(1);
        }
        node = node -> leftchild;
        (node -> parent) -> leftchild = NULL;
    }
    (*count)--;
    
    // place the last node on the top
    node -> parent = NULL;
    if(  node != rightchild ) {
        node -> rightchild = rightchild;
        if( (node -> rightchild) != NULL ) (node -> rightchild) -> parent = node;
    }
    else node -> rightchild = NULL;
    if(  node != leftchild ) {
        node -> leftchild = leftchild;
        if( (node -> leftchild) != NULL ) (node -> leftchild) -> parent = node;
    }
    else node -> leftchild = NULL;
    
    
    // find the right position for the last node in the tree after putting it on the top
    uleft = ( (node -> leftchild) != NULL ) ? (arc + ((node -> leftchild) -> ind)) -> cost : INFTY;
    uright = ( (node -> rightchild) != NULL ) ? (arc + ((node -> rightchild) -> ind)) -> cost : INFTY;
    u = (arc + (node -> ind )) -> cost;
    while( u > min(uleft, uright) ) {
        if( uleft < uright ) { // swap node with its left child
            swap(node,node -> leftchild,'l');
        }
        else if( uleft >= uright ) { // swap node and its rightchild
            swap(node,node -> rightchild,'r');
        }
        else {
            printf("Error in deltree\n");
            exit(1);
        }
        uleft = ( (node -> leftchild) != NULL ) ? (arc + ((node -> leftchild) -> ind)) -> cost : INFTY;
        uright = ( (node -> rightchild) != NULL ) ? (arc + ((node -> rightchild) -> ind)) -> cost : INFTY;
    }
    
    // find the new root
    while( node -> parent != NULL ) {
        node = node -> parent;
    }
    
    ind = node -> ind;
    
    return ind;
}


/***********   s w a p    a    n o d e     w i t h     i t s    c h i l d   ***********/
void swap( struct my_treenode *parent, struct my_treenode *child, char ch ) {
    // swap child and parent in the binary tree
    struct my_treenode *temp_leftchild, *temp_rightchild;
    
    
    temp_leftchild = child -> leftchild;
    temp_rightchild = child -> rightchild;
    child -> parent = parent -> parent;
    if( ch == 'l' ) { // child == left child
        child -> leftchild = parent;
        child -> rightchild = parent -> rightchild;
        if( (parent -> rightchild) != NULL ) (parent -> rightchild) -> parent = child;
    }
    else if( ch == 'r' ) { // child == right child
        child -> rightchild = parent;
        child -> leftchild = parent -> leftchild;
        if( (parent -> leftchild) != NULL )  (parent -> leftchild) -> parent = child;
    }
    if( child -> parent != NULL ) {
        if( parent == (child -> parent) -> leftchild ) (child -> parent) -> leftchild = child;
        else if(parent == (child -> parent) -> rightchild ) (child -> parent) -> rightchild = child;
    }
    
    parent -> parent = child;
    parent -> leftchild = temp_leftchild;
    parent -> rightchild = temp_rightchild;
    if( temp_leftchild != NULL ) (parent -> leftchild) -> parent = parent;
    if( temp_rightchild != NULL ) (parent -> rightchild) -> parent = parent;
    
    
}

/********************   u p d a t e      t r e e   *********************/
int updatetree( struct my_treenode *tree, struct my_treenode *node ) {
    double uleft, uright, u, uparent;
    int ind;
    
    uleft = ( (node -> leftchild) != NULL ) ? (arc + ((node -> leftchild) -> ind)) -> cost : INFTY;
    uright = ( (node -> rightchild) != NULL ) ? (arc + ((node -> rightchild) -> ind)) -> cost : INFTY;
    uparent = ( (node -> parent) != NULL ) ? (arc + ((node -> parent) -> ind)) -> cost : 0.0;
    u = (arc + (node -> ind)) -> cost;
    while( u > min(uleft, uright) ) {
        if( uleft < uright ) { // swap node with its left child
            swap(node,node -> leftchild,'l');
        }
        else if( uleft >= uright ) { // swap node and its rightchild
            swap(node,node -> rightchild,'r');
        }
        else {
            printf("Error in updatetree\n");
            exit(1);
        }
        uleft = ( (node -> leftchild) != NULL ) ? (arc + ((node -> leftchild) -> ind)) -> cost : INFTY;
        uright = ( (node -> rightchild) != NULL ) ? (arc + ((node -> rightchild) -> ind)) -> cost : INFTY;
    }
    while( node -> parent != NULL &&
          ((arc + (node -> ind)) -> cost) < ((arc + ((node -> parent) -> ind)) -> cost) ) {
        if( node == (node -> parent) -> leftchild ) {
            swap(node -> parent, node, 'l');
        }
        else if( node == (node -> parent) -> rightchild ) {
            swap(node -> parent, node, 'r');
        }
        else {
            printf("Error in updatetree\n");
            exit(1);
        }
    }
    // find the new root
    while( node -> parent != NULL ) node = node -> parent;
    tree = node;
    ind = tree -> ind;
    
    return ind;
}



/***************************************************/
void update_arcs(struct my_treenode *tree,double delta_u) {
    
    if(tree == NULL) {
        printf("tree = NULL\n");
    }
    else {
        ((arc + (tree -> ind)) -> cost) += delta_u ;
        if( tree -> leftchild != NULL ) update_arcs(tree -> leftchild,delta_u);
        if( tree -> rightchild != NULL ) update_arcs(tree -> rightchild,delta_u);
    }
}

/****************   m a k e      l i n k   *******************/

void make_link(int ind, char ch) {
    struct mylinklist *link;
    
    link = (struct mylinklist *)malloc(sizeof(struct mylinklist));
    
    link -> index = ind;
    link -> first = link;
    link -> last = link;
    link -> count = 1;
    link -> next = NULL;
    link -> prev = NULL;
    if( ch == 'l' ) {
        first[(arc +(link -> index)) -> head] = link;
        first[(arc +(link -> index)) -> tail] = link;
    }
    else if( ch == 'c' ) {
        cfirst[(arc +(link -> index)) -> head] = link;
        cfirst[(arc +(link -> index)) -> tail] = link;
    }
    else {
        printf("Error in make_link: ch = %c\n",ch);
        exit(1);
    }
    
}


/****************   m e r g e     l i n k l i s t s   ******************/

void merge_linklists(struct mylinklist *list1, struct mylinklist *list2, int ind, char ch) {
    // char == 'l': merge linklists of arcs
    // char == 'c': merge linklists of cycles
    struct mylinklist *link, *aux;
    
    link = (struct mylinklist *)malloc(sizeof(struct mylinklist));
    link -> index = ind;

    if( list1 != NULL  && list2 != NULL ) {
        if( (list1 -> count) < (list2 -> count) ) {
            aux = list1;
            list1 = list2;
            list2 = aux;
        }
        link -> first = list1;
        link -> prev = list1 -> last;
        (list1 -> last) -> next = link;
        link -> next = list2 -> first;
        link -> last = list2 -> last;
        list1 -> last = list2 -> last;
        list2  -> prev = link;
        list2 -> first = list1 -> first;
        if( ch == 'l' ) {
            do {
                first[(arc +(link -> index)) -> head] = list1;
                first[(arc +(link -> index)) -> tail] = list1;
                link -> first = list1;
                link = link -> next;
            }
            while( link != NULL );
        }
        else if( ch == 'c' ) {
            do {
                cfirst[(arc +(link -> index)) -> head] = list1;
                cfirst[(arc +(link -> index)) -> tail] = list1;
                link -> first = list1;
                link = link -> next;
            }
            while( link != NULL );
        }
        else {
            printf("Error in merge_linklists: ch = %c\n",ch);
            exit(1);
        }
        (list1 -> count) = (list1 -> count) + (list2 -> count) + 1;
    }
    else if( list1 != NULL && list2 == NULL ) {
        link -> first = list1;
        link -> prev = list1 -> last;
        (list1 -> last) -> next = link;
        list1 -> last = link;
        link -> next = NULL;
        (list1 -> count)++;
        if( ch == 'l' ) {
            first[(arc +(link -> index)) -> head] = list1;
            first[(arc +(link -> index)) -> tail] = list1;
        }
        else if( ch == 'c' ) {
            cfirst[(arc +(link -> index)) -> head] = list1;
            cfirst[(arc +(link -> index)) -> tail] = list1;
        }
        else {
            printf("Error in merge_linklists: ch = %c\n",ch);
            exit(1);
        }
    }
    else if( list1 == NULL && list2 != NULL ) {
        link -> first = list2;
        link -> prev = list2 -> last;
        (list2 -> last) -> next = link;
        list2 -> last = link;
        link -> next = NULL;
        (list2 -> count)++;
        if( ch == 'l' ) {
            first[(arc +(link -> index)) -> head] = list2;
            first[(arc +(link -> index)) -> tail] = list2;
        }
        else if( ch == 'c' ) {
            cfirst[(arc +(link -> index)) -> head] = list2;
            cfirst[(arc +(link -> index)) -> tail] = list2;
        }
        else {
            printf("Error in merge_linklists: ch = %c\n",ch);
            exit(1);
        }
    }
    else if( list1 == NULL && list2 == NULL ) {
        link -> first = link;
        link -> prev = NULL;
        link -> last = link;
        link -> next = NULL;
        link -> count = 1;
        if( ch == 'l' ) {
            first[(arc +(link -> index)) -> head] = link;
            first[(arc +(link -> index)) -> tail] = link;
        }
        else if( ch == 'c' ) {
            cfirst[(arc +(link -> index)) -> head] = link;
            cfirst[(arc +(link -> index)) -> tail] = link;
        }
        else {
            printf("Error in merge_linklists: ch = %c\n",ch);
            exit(1);
        }
    }
    
}



/****************   f i n d     c y c l e   *****************/

int find_cycle(int tail, int head, int newarc) {
    int  root, tail0, i, icycle;
    double delta_u;
    
    tail0 = tail;
    root = next_arc_index[tail0];
    
    while( min_arc_index[head] != newarc ) {
        // merge linklists of cycles and the corresponding trees of unused outgoing arcs
        merge_linklists(cfirst[tail0],cfirst[head],min_arc_index[tail],'c');
        delta_u = ((arc + newarc) -> cost) - ((arc + min_arc_index[head]) -> cost);
        
        if( next_arc_index[head] >= 0 ) { /* if the tree of unused outgoing arcs is nonempty */
            update_arcs(treenode + next_arc_index[head],delta_u);
            if( root >= 0 ) { // tail has nonempty tree of outgoing arcs
                // add tree with fewer nodes to the one with more nodes
                if( count_vertex[tail0] >= count_vertex[head] ) {
                    root = addtree(treenode + root,treenode + next_arc_index[head],count_vertex + tail0);
                    count_vertex[head] = count_vertex[tail0];
                }
                else {
                    root = addtree(treenode + next_arc_index[head],treenode + root,count_vertex + head);
                    count_vertex[tail0] = count_vertex[head];
                }
            }
            else { // tail has no outgoing arcs left in the set B_{tail}
                root = next_arc_index[head];
                count_vertex[tail0] = count_vertex[head];
            }
        }
        tail = head;
        head = (arc + min_arc_index[head]) -> head;
    } /* end while */
    

    return root;
    
}

/******************    m e r g e     v l i s t s   ******************/

void merge_vlists(int head,int tail) {
    struct simplelist *vlist;
    
    (vfirst[head] -> last) -> next = vfirst[tail];
    vfirst[head] -> last = vfirst[tail] -> last;
    vlist = vfirst[tail];
    while( vlist != NULL ) {
        vfirst[vlist -> index] = vfirst[head];
        vlist = vlist -> next;
    }
}


/*****************   f i n d      o p t i m a l     w - g r a p h   *******************/

void find_optimal_wgraph(int sink) {
    struct simplelist *incoming;

    incoming = infirst[sink];
    used[sink] = 'y';
    while( incoming != NULL ) {
        if( used[(arc + (incoming -> index)) -> tail ] == 'n') {
            if( WRITE_FILE == 'y' ) {
                fprintf(fg,"%i\t%i\n",(arc + (incoming -> index)) -> tail,(arc + (incoming -> index)) -> head);
            }
            cc_count++;
            find_optimal_wgraph((arc + (incoming -> index)) -> tail);
        }
        incoming = incoming -> next;
    }
}

/*****************   m a x i m a l     c o n n e c t e d      c o m p o n e n t s   *******************/
void maximal_connected_components(void) {
    // FIND THE MAXIMAL CONNECTED COMPONENTS
    int isinkmax, nmax, isink, k;
    struct simplelist *vlist;
    char fname[100];
    FILE *fss;
    
    printf("In maximal_connected_components()\n");
    printf("In total, %i connected components\n",n_connected_components);
    isinkmax = 0;
    nmax = 0;
    if( n_connected_components > 1 ) {
        fss = fopen("sinks_and_sizes.txt","w");
    }
    for( isink = 0; isink < n_connected_components; isink++ ) {
        if( WRITE_FILE == 'y' ) {
            sprintf(fname,"owg%i.txt",isink);
            fg = fopen(fname,"w");
        }
        for( k = 0; k < nstates; k++ ) used[k] = 'n';
        cc_count = 1;
        find_optimal_wgraph(vfirst[main_sink[isink]] -> index);
        if( WRITE_FILE == 'y' ) fclose(fg);
        size_con_comp[isink] = cc_count;
        printf("sink[%i] = %i, size = %i\n",isink,main_sink[isink],cc_count);
        if( cc_count > nmax ) {
            nmax = cc_count;
            isinkmax = isink;
        }
        if( n_connected_components > 1 ) {
            fprintf(fss,"%i\t%i\n",main_sink[isink],cc_count);
        }
    }
    if( n_connected_components > 1 ) {
        fclose(fss);
        printf("The maximal connected component: sink = %i,\t size = %i\n",main_sink[isinkmax],size_con_comp[isinkmax]);
        fss = fopen("max_con_comp.txt","w");
        vlist = vfirst[main_sink[isinkmax]];
        while( vlist != NULL ) {
            fprintf(fss,"%i\n",vlist -> index );
            vlist = vlist -> next;
        }
        fclose(fss);
    }
    
}



/**************** M A I N *****************/

int main() {
    clock_t CPUbegin;
    double cpu;
    char ch;
    long npath;
    
    readinput();
    CPUbegin = clock();
    init();
    
    add_arc();
    
    cpu=(clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
    printf("cputime = %g\n",cpu);
    
    return 0;
}





