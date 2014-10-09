#include <R.h>
#include <Rinternals.h>

static int iii;


void net_reorder(int node, int n, int *e1, int *e2, int *neworder, int *L, int *pos1, int *pos2, int *visited)
{   
    int newnode=0;
//    node = node;
    int i = pos1[node], j, k;
    visited[node]=1L;
    for (j = 0; j < pos2[node]; j++) {
        k = L[i + j];
		neworder[iii++] = k + 1;
        newnode = e2[k];
		if ((newnode > n) & (visited[newnode] == 0L))/* is it an internal edge or visited before */
			net_reorder(e2[k], n, e1, e2, neworder, L, pos1, pos2, visited);
	}
}


void order_networx(int *ntips, int *nEdges, int *mNodes, int *e1, int *e2, int *root, int * neworder){
    int i, j, m = mNodes[0]+1L;//m=*mNodes - *ntips;
    int *pos1, *pos2, *tmp, *L; 
    tmp = (int*)R_alloc(m, sizeof(int));
    pos1 = (int*)R_alloc(m, sizeof(int));
    pos2 = (int*)R_alloc(m, sizeof(int));
    L = (int*)R_alloc(*nEdges, sizeof(int));
    for(i = 0; i<m; i++) tmp[i]=0; 
    for(i = 0; i<m; i++) pos1[i]=0; 
    for(i = 0; i<m; i++) pos2[i]=0;

    for(i = 0; i<*nEdges; i++) tmp[e1[i]]++;
  // pos1[0]=0;
    for(i = 0; i<(*mNodes); i++) pos1[i+1]=pos1[i] + tmp[i];
     for(i = 0; i<*nEdges; i++){
         j = e1[i];
         L[pos1[j] + pos2[j]] = i;  
         pos2[j]++;
    }
    for(i = 0; i<m; i++) pos2[i]=0L;
    net_reorder(*root, *ntips, e1, e2, neworder, L, pos1, tmp, pos2);
//   for(i = 0; i<*nEdges; i++) neworder[i] = L[i];
}




