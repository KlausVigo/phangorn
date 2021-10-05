/*
 * sprdist.c
 *
 * (c) 2016-2019 Leonardo de Oliveira Martins (leomrtns@gmail.com)
 *
 *
 * This code may be distributed under the GNU GPL
 *
 */

#include <Rmath.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>     /* standard integer types (int32_t typedef etc.) [C99]*/


#define true  1U /*!< Boolean TRUE  */
#define false 0U /*!< Boolean FALSE */

typedef unsigned char bool;
typedef struct splitset_struct* splitset;
typedef struct hungarian_struct* hungarian;
typedef struct bipartition_struct* bipartition;
typedef struct bipsize_struct* bipsize;

struct splitset_struct
{
  int size, spsize, spr, spr_extra, rf, hdist; /*! \brief spr, extra prunes for spr, rf distances and hdist=assignment cost */
  int n_g, n_s, n_agree, n_disagree;
  bipartition *g_split, *s_split, *agree, *disagree;
  bipartition prune;
  hungarian h; /* hungarian method for solving the assignment between edges */
  bool match;  /*! \brief do we want to calculate the minimum cost assignment */
};

struct hungarian_struct
{
  int **cost, *col_mate; /*! \brief cost matrix, and col_mate[row] with column match for row */
  int size,  /*! \brief assignment size. Cost is a square matrix, so size should be an overestimate where "missing" nodes are added w/ cost zero */
      initial_cost, /*! \brief sum of lowest input cost values for each column. The hungarian method rescales them so that minimum per column is zero */
      final_cost;   /*! \brief our final cost is on rescaled cost matrix, therefore to restore the "classical" optimal cost one should sum it with initial_cost */
  int *unchosen_row, *row_dec, *slack_row, *row_mate, *parent_row, *col_inc, *slack; /* aux vectors */
};

/*! \brief Bit-string representation of splits. */
struct bipartition_struct
{
  unsigned long long *bs;  /*! \brief Representation of a bipartition by a vector of integers (bitstrings). */
  int n_ones;  /*! \brief Counter (number of "one"s) */
  bipsize n;  /*! \brief number of bits (leaves), vector size and mask */
  int ref_counter;  /*! \brief How many times this struct is being referenced */
};

struct bipsize_struct
{
  unsigned long long mask;/*! \brief mask to make sure we consider only active positions (of last bitstring) */
  int ints, bits, original_size;  /*! \brief Vector size and total number of elements (n_ints = n_bits/(8*sizeof(long long)) +1). */
  int ref_counter;  /*! \brief How many times this struct is being referenced */
};


/*! \brief Allocate space for splitset structure (two vectors of bipartitions), for simple comparisons */
splitset new_splitset (int nleaves, int nsplits);
/*! \brief free memory allocated for splitset structure */
void del_splitset (splitset split);
/*! \brief low level function that does the actual SPR and hdist calculation based on a filled splitset struct */
int dSPR_topology_lowlevel (splitset split);
/*! \brief function used by qsort for a vector of bipartitions (from smaller to larger) */
int compare_splitset_bipartition_increasing (const void *a1, const void *a2);

/* BELOW: low level functions that work with bipartitions */
void split_create_agreement_list (splitset split);
void split_remove_agree_edges (splitset split, bipartition *b, int *nb);
void split_remove_duplicates (bipartition *b, int *nb);
void split_compress_agreement (splitset split);
void split_create_disagreement_list (splitset split);
void split_disagreement_assign_match (splitset split);
void split_find_small_disagreement (splitset split);
void split_remove_small_disagreement (splitset split);
void split_minimize_subtrees (splitset split);
void split_remove_redundant_bit (splitset split, int id);
void split_replace_bit (splitset split, int to, int from);
void split_new_size (splitset split, int size, bool update_bipartitions);
void split_swap_position (bipartition *b, int i1, int i2);

/* BELOW: Hungarian method for bipartite matching (assignment) */
hungarian new_hungarian (int size);
void hungarian_reset (hungarian p);
void hungarian_update_cost (hungarian p, int row, int col, int cost);
void del_hungarian (hungarian p);
void hungarian_solve (hungarian p, int this_size);

/* BELOW: Memory-efficient, fast bipartition comparisions based on 64bit representation of splits */
bipartition new_bipartition (int size);
bipsize new_bipsize (int size);
bipartition new_bipartition_copy_from (const bipartition from);
bipartition new_bipartition_from_bipsize (bipsize n);
void del_bipartition (bipartition bip);
void del_bipsize (bipsize n);
void bipsize_resize (bipsize n, int nbits);
void bipartition_initialize (bipartition bip, int position);
void bipartition_zero (bipartition bip);
void bipartition_set (bipartition bip, int position);
void bipartition_set_lowlevel (bipartition bip, int i, int j);
void bipartition_unset (bipartition bip, int position);
void bipartition_unset_lowlevel (bipartition bip, int i, int j);
void bipartition_copy (bipartition to, const bipartition from);
void bipartition_OR (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
void bipartition_AND (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
void bipartition_ANDNOT (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
void bipartition_XOR (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
void bipartition_XORNOT (bipartition result, const bipartition b1, const bipartition b2, bool update_count);
void bipartition_NOT (bipartition result, const bipartition bip);
void bipartition_count_n_ones (const bipartition bip);
void bipartition_to_int_vector (const bipartition b, int *id, int vecsize);
bool bipartition_is_equal (const bipartition b1, const bipartition b2);
bool bipartition_is_equal_bothsides (const bipartition b1, const bipartition b2);
bool bipartition_is_larger (const bipartition b1, const bipartition b2);
void bipartition_flip_to_smaller_set (bipartition bip);
bool bipartition_is_bit_set (const bipartition bip, int position);
bool bipartition_contains_bits (const bipartition b1, const bipartition b2);
// void bipartition_print_to_stdout (const bipartition b1);
void bipartition_replace_bit_in_vector (bipartition *bvec, int n_b, int to, int from, bool reduce);
void bipartition_resize_vector (bipartition *bvec, int n_b);

/*! \brief Main SPR calculation function, to be used within R */
SEXP C_sprdist (SEXP bp1, SEXP bp2, SEXP lt) {
  int i, j, n_leaves = INTEGER(lt)[0];
  SEXP result;
  double *res;
  splitset split;

  PROTECT(result = allocVector(REALSXP, 4));
  res = REAL(result);

//  if (length(bp1) != length(bp2)) error ("number of bipartitions given to C_sprdist are not the same");

  split = new_splitset (n_leaves, length(bp1));
  for (i=0; i < split->size; i++) {
//    for (j=0; j < length(VECTOR_ELT (bp1, i)); j++) printf (";;%d  ", INTEGER (VECTOR_ELT (bp1, i))[j]);
    for (j=0; j < length(VECTOR_ELT (bp1, i)); j++) bipartition_set (split->g_split[i], INTEGER (VECTOR_ELT (bp1, i))[j] - 1);
    for (j=0; j < length(VECTOR_ELT (bp2, i)); j++) bipartition_set (split->s_split[i], INTEGER (VECTOR_ELT (bp2, i))[j] - 1);
  }
  dSPR_topology_lowlevel (split);
  res[0] = split->spr;
  res[1] = split->spr_extra;
  res[2] = split->rf;
  res[3] = split->hdist;
  del_splitset (split);

  UNPROTECT(1); // result
  return(result);
}

/* functions below should not be called outside this scope */

splitset
new_splitset (int nleaves, int nsplits)
{
  splitset split;
  int i;

  split = (splitset) malloc (sizeof (struct splitset_struct));
  split->n_g = split->n_s = split->size = nsplits;
  split->n_agree = split->n_disagree = 0;
  split->prune = NULL;
  split->match = true; /* do we want to calculate the assignment matching cost (using hungarian() )? */
  split->spr = split->spr_extra = split->rf = split->hdist = 0;

  split->g_split = (bipartition*) malloc (split->size * sizeof (bipartition));
  split->s_split = (bipartition*) malloc (split->size * sizeof (bipartition));
  split->g_split[0] = new_bipartition (nleaves);
  split->s_split[0] = new_bipartition (nleaves);
  for (i = 1; i < split->size; i++) {
    split->g_split[i] = new_bipartition_from_bipsize (split->g_split[0]->n); /* use same bipsize */
    split->s_split[i] = new_bipartition_from_bipsize (split->s_split[0]->n);
  }


  split->agree    = (bipartition*) malloc (split->size * sizeof (bipartition));
  split->disagree = (bipartition*) malloc (split->size * split->size * sizeof (bipartition));
  split->agree[0]    = new_bipartition (nleaves); // this bipsize will be recycled below
  split->disagree[0] = new_bipartition (nleaves);
  for (i = 1; i < split->size; i++)               split->agree[i]    = new_bipartition_from_bipsize (split->agree[0]->n);
  for (i = 1; i < split->size * split->size; i++) split->disagree[i] = new_bipartition_from_bipsize (split->disagree[0]->n);
  split->prune = new_bipartition_from_bipsize (split->disagree[0]->n);

  split->h = new_hungarian (split->size);

  return split;
}

void
del_splitset (splitset split)
{
  int i;
  if (!split) return;

  del_bipartition (split->prune);
  if (split->disagree) {
    for (i = split->size * split->size - 1; i >= 0; i--) del_bipartition (split->disagree[i]);
    free (split->disagree);
  }
  if (split->agree) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->agree[i]);
    free (split->agree);
  }
  if (split->g_split) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->g_split[i]);
    free (split->g_split);
  }
  if (split->s_split) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->s_split[i]);
    free (split->s_split);
  }
  del_hungarian (split->h);
  free (split);
}

int
compare_splitset_bipartition_increasing (const void *a1, const void *a2)
{ /* similar to bipartition_is_larger() */
  bipartition *b1 = (bipartition *) a1;
  bipartition *b2 = (bipartition *) a2;
  int i;

  if ((*b1)->n_ones > (*b2)->n_ones) return 1;
  if ((*b1)->n_ones < (*b2)->n_ones) return -1;

  for (i = (*b1)->n->ints - 1; (i >= 0) && ((*b1)->bs[i] == (*b2)->bs[i]); i--); /* find position of distinct bipartition elem*/
  if (i < 0) return 0; /* identical bipartitions */
  if ((*b1)->bs[i] > (*b2)->bs[i]) return 1;
  else return -1;
}


int
dSPR_topology_lowlevel (splitset split)
{
  int i = 0, mismatch = -1;

  for (i=0; i < split->size; i++) {
    bipartition_flip_to_smaller_set (split->g_split[i]);
    bipartition_flip_to_smaller_set (split->s_split[i]);
  }
  qsort (split->g_split, split->size, sizeof (bipartition), compare_splitset_bipartition_increasing);
  qsort (split->s_split, split->size, sizeof (bipartition), compare_splitset_bipartition_increasing);
  //for (i = 0; i < split->n_g; i++) bipartition_print_to_stdout (split->g_split[i]); printf ("G   ::DEBUG 0 ::\n");
  //for (i = 0; i < split->n_s; i++) bipartition_print_to_stdout (split->s_split[i]); printf ("S\n");

  i++; /* to trick -Werror, since we don't use it unless for debug */
  while (mismatch) {
    split_create_agreement_list (split);  // vector of identical bipartitions
    split_compress_agreement (split);     // iterative replacement of cherry by new leaf

//    for (i = 0; i < split->n_g; i++) bipartition_print_to_stdout (split->g_split[i]); printf ("G   ::DEBUG::\n");
//    for (i = 0; i < split->n_s; i++) bipartition_print_to_stdout (split->s_split[i]); printf ("S\n");
//    for (i = 0; i < split->n_agree; i++) bipartition_print_to_stdout (split->agree[i]); printf ("A\n");

    if (mismatch == -1) split->rf = split->n_g + split->n_s;
    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
    if (!mismatch) return split->spr;

    split_create_disagreement_list (split); // vector of smallest disagreements
    split_disagreement_assign_match (split); /* assignment matching between edges using hungarian method (split->hdist after first time) */

    split_remove_duplicates (split->disagree, &(split->n_disagree)); // some elements are equal; this function also qsorts
    split_find_small_disagreement (split);  // could also be one leaf only

    //for (i = 0; i < split->n_disagree; i++) { bipartition_print_to_stdout (split->disagree[i]); printf ("\n"); }
    //printf ("{%d} prune: ", split->n_disagree); bipartition_print_to_stdout (split->prune); printf ("\n");

    split->spr++;
    split_remove_small_disagreement (split);

    split_minimize_subtrees (split);
    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
  }
  return split->spr;
}

void
split_create_agreement_list (splitset split)
{
  int s, g;
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)
    if (bipartition_is_equal (split->g_split[g], split->s_split[s])) {
      bipartition_copy (split->agree[split->n_agree++], split->g_split[g]);
      split->n_g--; split_swap_position (split->g_split, g, split->n_g); /* if we don't swap them, we lose ref to "old" value on g_split[] */
      split->n_s--; split_swap_position (split->s_split, s, split->n_s);
      g--; s = split->n_s; /* pretend loop finished, examine again with new values */
    }
  split_remove_agree_edges (split, split->g_split, &(split->n_g));
  split_remove_agree_edges (split, split->s_split, &(split->n_s));
}

void
split_remove_agree_edges (splitset split, bipartition *b, int *nb)
{
  int i, a;
  for (i = 0; i < (*nb); i++) for (a = 0; a < split->n_agree; a++)
    if (bipartition_is_equal (b[i], split->agree[a])) {
      (*nb)--;
      split_swap_position (b, i, (*nb));
      i--;
      a = split->n_agree; /* loop again over new value */
    }
}

void
split_remove_duplicates (bipartition *b, int *nb)
{
  int i, j;
  bipartition pivot;
  if ((*nb) < 2) return; /* only if we have a vector with > 1 element */
  qsort (b, (*nb), sizeof (bipartition), compare_splitset_bipartition_increasing);

  for (i = (*nb) - 1; i >= 1; i--)
    if (bipartition_is_equal (b[i], b[i-1])) {
      pivot = b[i]; /* do not lose a pointer to this element */
      for (j = i; j < (*nb)-1; j++) b[j] = b[j+1];
      b[j] = pivot; /* j = (*nb) - 1, which will become obsolete through next line -->  (*nb)-- */
      (*nb)--;
    }
}

void
split_compress_agreement (splitset split)
{
  int i, j, pair[2];

  for (i = 0; i < split->n_agree; i++) if (split->agree[i]->n_ones == 2) { /* cherry in common, can be represented by just one leaf */
    bipartition_to_int_vector (split->agree[i], pair, 2);
    split_remove_redundant_bit (split, pair[1]);
    split_new_size (split,split->agree[0]->n->bits - 1, false); /* false = do not recalculate every bipartition's last elem */
    bipartition_resize_vector (split->agree, split->n_agree);
    for (j = 0; j < split->n_agree; j++) { /* minimize subtree size and remove single leaves */
      bipartition_flip_to_smaller_set (split->agree[j]); /* agree only */
      if (split->agree[j]->n_ones < 2) split_swap_position (split->agree, j--, --split->n_agree);
    }
    i = -1; /* redo all iterations, with new info (agree[] will be smaller) */
  }
  bipartition_resize_vector (split->g_split, split->n_g);
  bipartition_resize_vector (split->s_split, split->n_s);
}

void
split_create_disagreement_list (splitset split)
{
  int g, s;

  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++) {
    bipartition_XOR (split->disagree[g * split->n_s + s], split->g_split[g], split->s_split[s], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[g * split->n_s + s]);
  }
  split->n_disagree = split->n_g * split->n_s;
}

void
split_disagreement_assign_match (splitset split)
{ /* also calculates split->hdist */
  int g, s, max_n, sum = 0;

  if (split->n_g > split->n_s) max_n = split->n_g;
  else                         max_n = split->n_s;
  if (max_n < 2) return;

  hungarian_reset (split->h);
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)
    hungarian_update_cost (split->h, g, s, split->disagree[g * split->n_s + s]->n_ones);
  hungarian_solve (split->h, max_n);
  /* now split->h->col_mate will have the pairs */
  /* if we do the matching below it becomes much faster, but we may miss the best prune subtrees in a few cases (do not compromise the algo) */
  split->n_disagree = 0;
  for (g = 0; g < max_n; g++) if ((g < split->n_g) && ( split->h->col_mate[g] < split->n_s)) { /* some matchings might be to dummy edges */
    bipartition_XOR (split->disagree[split->n_disagree], split->g_split[g], split->s_split[split->h->col_mate[g]], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[split->n_disagree++]);
    sum += split->disagree[split->n_disagree-1]->n_ones;
  }
  if (split->match) { split->hdist = split->h->final_cost+split->h->initial_cost; split->match = false; }
}

void
split_find_small_disagreement (splitset split)
{
  bipartition dis;
  int a, d;

  bipartition_copy (split->prune, split->disagree[0]); /* smallest, in case we don't find a better one in loop below */
  if (split->prune->n_ones < 2) return;

  dis = new_bipartition_from_bipsize (split->disagree[0]->n);
  for (d = 0; d < split->n_disagree; d++) for (a = 0; a < split->n_agree; a++) {
    if ((split->disagree[d]->n_ones == split->agree[a]->n_ones) ||
        (split->disagree[d]->n_ones == (split->agree[a]->n->bits - split->agree[a]->n_ones))) {
      bipartition_XOR (dis, split->disagree[d], split->agree[a], true);
      if      (!dis->n_ones)               { bipartition_copy (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
      else if (dis->n_ones == dis->n->bits) { bipartition_NOT (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
    }
  }
  /* check if prune nodes are all on same side of a tree or if they are actually two SPRs (one from each tree) */
  for (d = 0; d < split->n_g; d++) {
    if (!bipartition_contains_bits (split->g_split[d], split->prune)) {
      bipartition_NOT (dis, split->g_split[d]);
      if (!bipartition_contains_bits (dis, split->prune)) { split->spr_extra++; d = split->n_g; }
    }
  }
  del_bipartition (dis);
}

void
split_remove_small_disagreement (splitset split)
{
  int *index, i, j = split->prune->n_ones - 1, k = 0, size = split->agree[0]->n->bits;

  index = (int*) malloc (split->prune->n_ones * sizeof (int));
  bipartition_to_int_vector (split->prune, index, split->prune->n_ones);

  for (i = size - 1; i >= (size - split->prune->n_ones); i--) {
    if (index[k] >= (size - split->prune->n_ones)) i = -1;
    else {
      if (i == index[j]) j--;
      else split_replace_bit (split, index[k++], i);
    }
  }

  split_new_size (split,size - split->prune->n_ones, true);
  if (index) free (index);
}

void
split_minimize_subtrees (splitset split)
{
  int i;

  for (i = 0; i < split->n_s; i++) {
    bipartition_flip_to_smaller_set (split->s_split[i]);
    if (split->s_split[i]->n_ones < 2) { split->n_s--; split_swap_position (split->s_split, i, split->n_s); i--; }
  }
  for (i = 0; i < split->n_g; i++) {
    bipartition_flip_to_smaller_set (split->g_split[i]);
    if (split->g_split[i]->n_ones < 2) { split->n_g--; split_swap_position (split->g_split, i, split->n_g); i--; }
  }
  for (i = 0; i < split->n_agree; i++) {
    bipartition_flip_to_smaller_set (split->agree[i]);
    if (split->agree[i]->n_ones < 2) { split->n_agree--; split_swap_position (split->agree, i, split->n_agree); i--; }
  }
}

void
split_remove_redundant_bit (splitset split, int id)
{
  int  last = split->agree[0]->n->bits-1;
  if (id < last) split_replace_bit (split, id, last);
}

void
split_replace_bit (splitset split, int to, int from)
{
  if (from <= to) return;
  /*not needed for disagree[] */
  bipartition_replace_bit_in_vector (split->agree,   split->n_agree, to, from, true);
  bipartition_replace_bit_in_vector (split->g_split, split->n_g,     to, from, true);
  bipartition_replace_bit_in_vector (split->s_split, split->n_s,     to, from, true);
}

void
split_new_size (splitset split, int size, bool update_bipartitions)
{
  bipsize_resize (split->g_split[0]->n, size);
  bipsize_resize (split->s_split[0]->n, size);
  bipsize_resize (split->agree[0]->n, size);
  bipsize_resize (split->disagree[0]->n, size);
  if (update_bipartitions) {
    bipartition_resize_vector (split->g_split, split->n_g);
    bipartition_resize_vector (split->s_split, split->n_s);
    bipartition_resize_vector (split->agree, split->n_agree);
  }
}

void
split_swap_position (bipartition *b, int i1, int i2)
{
  bipartition pivot = b[i1];
  b[i1] = b[i2];
  b[i2] = pivot;
}


/* The hungarian method below is copied from http://www.informatik.uni-freiburg.de/~stachnis/misc.html
 * The (edited) original message follows:
 *
 ** libhungarian by Cyrill Stachniss, 2004  Solving the Minimum Assignment Problem using the
 ** Hungarian Method.         ** This file may be freely copied and distributed! **
 **
 ** Parts of the used code was originally provided by the "Stanford GraphGase", but I made changes to this code.
 ** As asked by  the copyright node of the "Stanford GraphGase", I hereby proclaim that this file are *NOT* part of the
 ** "Stanford GraphGase" distrubition! */

void
hungarian_reset (hungarian p)
{
  int i, j;

  for (i = 0; i < p->size; i++) {
    p->col_mate[i] = p->unchosen_row[i] = p->row_dec[i] = p->slack_row[i] = p->row_mate[i] = p->parent_row[i] = p->col_inc[i] = p->slack[i] = 0;
    for (j = 0; j < p->size; j++) p->cost[i][j] = 0;
  }
  p->final_cost = 0;
}

hungarian
new_hungarian (int size)
{
  int i;
  hungarian p;

  p = (hungarian) malloc (sizeof (struct hungarian_struct));
  p->size = size; /* n_rows = n_columns; if it's not, fill with zeroes (no cost) */
  p->cost = (int**) malloc (size * sizeof (int*));
  for (i = 0; i < p->size; i++)
    p->cost[i] = (int*) malloc (size * sizeof (int));
  /* edges would be assignment_matrix[ i * ncols + col_mate[i] ] = true; and other elems "false" (but we don't use the matrix notation) */
  p->col_mate     = (int*) malloc (size * sizeof (int)); /* for a given row node, col_mate[row] is the assigned col node */
  p->unchosen_row = (int*) malloc (size * sizeof (int));
  p->row_dec      = (int*) malloc (size * sizeof (int));
  p->slack_row    = (int*) malloc (size * sizeof (int));
  p->row_mate     = (int*) malloc (size * sizeof (int));
  p->parent_row   = (int*) malloc (size * sizeof (int));
  p->col_inc      = (int*) malloc (size * sizeof (int));
  p->slack        = (int*) malloc (size * sizeof (int));

  hungarian_reset (p);
  return p;
}

void
hungarian_update_cost (hungarian p, int row, int col, int cost)
{
  if (row >= p->size) return;
  if (col >= p->size) return;
  p->cost[row][col] = cost;
}

void
del_hungarian (hungarian p)
{
  int i;
  if (!p) return;
  if (p->cost) {
    for (i = p->size - 1; i >= 0; i--) if (p->cost[i]) free (p->cost[i]);
    free (p->cost);
  }
  free (p->col_mate); /* this is the important one, with i assigned to col_mate[i] */
  free (p->slack);
  free (p->col_inc);
  free (p->parent_row);
  free (p->row_mate);
  free (p->slack_row);
  free (p->row_dec);
  free (p->unchosen_row);
  free (p);
}

void
hungarian_solve (hungarian p, int this_size)
{
  int i, j, nrows = this_size, ncols = this_size, k, l, s, t, q, unmatched;
  p->final_cost = p->initial_cost = 0;

  if (this_size > p->size) { p->final_cost = -1; return; } /* we don't call biomcmc_error(), but it *is* an error! */

  for (l = 0; l < ncols; l++) { // Begin subtract column minima in order to start with lots of zeroes 12
    s = p->cost[0][l];
    for (k = 1; k < nrows; k++) if (p->cost[k][l] < s) s = p->cost[k][l];
    p->initial_cost += s; /* this should be added to final_cost to have classical assignment cost; here we distinguish them */
    if (s!=0)	for (k = 0; k < nrows; k++) p->cost[k][l] -= s;
  } // End subtract column minima in order to start with lots of zeroes 12

  // Begin initial state 16
  t=0;
  for (l = 0; l < ncols; l++)  { // n => num_cols
    p->row_mate[l]= -1;
    p->parent_row[l]= -1;
    p->col_inc[l]=0;
    p->slack[l]= 0x7FFFFFFF;
  }
  for (k = 0; k < nrows; k++) { // m => num_rows
    s = p->cost[k][0];
    for (l = 1; l < ncols; l++) if (p->cost[k][l] < s) s = p->cost[k][l];
    p->row_dec[k]=s;
    for (l = 0; l < ncols; l++) if ((s==p->cost[k][l]) && (p->row_mate[l] < 0)) {
      p->col_mate[k] = l;
      p->row_mate[l] = k;  // fprintf(stderr, "matching col %d==row %d\n",l,k);
      goto row_done;
    }
    p->col_mate[k] = -1;  // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    p->unchosen_row[t++] = k;
row_done:
    ;
  }
  // End initial state 16

  // Begin Hungarian algorithm 18
  if (t==0)    goto done;
  unmatched=t;
  while (1) {
    q=0; // fprintf(stderr, "Matched %d rows.\n",m-t);
    while (1)	{
      while (q<t) {
         { // Begin explore node q of the forest 19
          k = p->unchosen_row[q];
          s=p->row_dec[k];
          for (l=0;l<ncols;l++) if (p->slack[l]) {
            int del;
            del = p->cost[k][l] - s + p->col_inc[l];
            if (del < p->slack[l]) {
              if (del==0) {
                if (p->row_mate[l]<0)  goto breakthru;
                p->slack[l]=0;
                p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n", t,row_mate[l],l,k);
                p->unchosen_row[t++]=p->row_mate[l];
              }
              else { p->slack[l]=del; p->slack_row[l]=k; }
            }
          }
         } // End explore node q of the forest 19
        q++;
      }

      // Begin introduce a new zero into the matrix 21
      s = 0x7FFFFFFF;
      for (l = 0;l < ncols; l++) if (p->slack[l] && p->slack[l] < s) s = p->slack[l];
      for (q = 0; q < t; q++) p->row_dec[ p->unchosen_row[q] ] += s;
      for (l = 0; l < ncols; l++) if (p->slack[l]) {
        p->slack[l]-=s;
        if (p->slack[l]==0) {  // Begin look at a new zero 22
          k = p->slack_row[l]; // fprintf(stderr, "Decreasing uncovered elements by %d produces zero at [%d,%d]\n", s,k,l);
          if (p->row_mate[l]<0)  {
            for (j=l+1;j<ncols;j++)  if (p->slack[j]==0) p->col_inc[j]+=s;
            goto breakthru;
          }
          else {
            p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k);
            p->unchosen_row[t++]=p->row_mate[l];
          }
        } // End look at a new zero 22

      }
      else  p->col_inc[l]+=s;
      // End introduce a new zero into the matrix 21
    }
breakthru:
    // fprintf(stderr, "Breakthrough at node %d of %d!\n",q,t);
    while (1)	{    // Begin update the matching 20
      j=p->col_mate[k];
      p->col_mate[k]=l;
      p->row_mate[l]=k; // fprintf(stderr, "rematching col %d==row %d\n",l,k);
      if (j<0)   break;
      k=p->parent_row[j];
      l=j;
    }    // End update the matching 20
    if (--unmatched==0)	goto done;
    // Begin get ready for another stage 17
    t=0;
    for (l=0;l<ncols;l++) {
      p->parent_row[l]= -1;
      p->slack[l]=0x7FFFFFFF;
    }
    for (k=0;k<nrows;k++) if (p->col_mate[k]<0) p->unchosen_row[t++]=k; // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    // End get ready for another stage 17
  }
done:
  // Begin doublecheck the solution 23
  for (k = 0; k < nrows; k++) for (l=0;l<ncols;l++) if (p->cost[k][l] < p->row_dec[k] - p->col_inc[l]) { p->final_cost = -1;  return;} //printf ("\n**\n");
  for (k = 0; k < nrows; k++) {
    l=p->col_mate[k];
    if ((l < 0) || (p->cost[k][l] != p->row_dec[k] - p->col_inc[l])) { p->final_cost = -1; return; }
  }
  k=0;
  for (l=0;l<ncols;l++) if (p->col_inc[l])  k++;
  if (k>nrows) { p->final_cost = -1; return; }
  // End doublecheck the solution 23
  // End Hungarian algorithm 18

  for (k = 0; k < nrows; ++k) for (l = 0; l < ncols; ++l) p->cost[k][l] = p->cost[k][l] - p->row_dec[k] + p->col_inc[l];
  for (i = 0; i < nrows; i++) p->final_cost += p->row_dec[i];
  for (i = 0; i < ncols; i++) p->final_cost -= p->col_inc[i]; // fprintf(stderr, "Cost is %d\n",cost);
}
 /* BELOW are the bipartition functions (memory-efficient storage of bipartitions on 64 bits */

int BitStringSize = 8 * sizeof (unsigned long long);

bipartition
new_bipartition (int size)
{
  bipartition bip;
  int i;

  bip = (bipartition) malloc (sizeof (struct bipartition_struct));
  bip->n = new_bipsize (size);
  bip->n_ones = 0;
  bip->ref_counter = 1;

  bip->bs = (unsigned long long*) malloc (bip->n->ints * sizeof (unsigned long long));
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = 0ULL;

  return bip;
}

bipsize
new_bipsize (int size)
{
  bipsize n;
  int i;

  n = (bipsize) malloc (sizeof (struct bipsize_struct));
  n->bits = n->original_size = size;
  n->ref_counter = 1;
  n->ints = size/BitStringSize + 1;
  n->mask = 0ULL;
  for (i=0; i < n->bits%BitStringSize; i++) n->mask |= (1ULL << i); /* disregard other bits */

  return n;
}

bipartition
new_bipartition_copy_from (const bipartition from)
{
  bipartition bip;
  int i;

  bip = (bipartition) malloc (sizeof (struct bipartition_struct));
  bip->n = new_bipsize (from->n->bits);
  bip->n_ones = from->n_ones;
  bip->ref_counter = 1;

  bip->bs = (unsigned long long*) malloc (bip->n->ints * sizeof (unsigned long long));
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = from->bs[i];

  return bip;
}

bipartition
new_bipartition_from_bipsize (bipsize n)
{
  bipartition bip;
  int i;

  bip = (bipartition) malloc (sizeof (struct bipartition_struct));
  bip->n = n;
  bip->n->ref_counter++;
  bip->n_ones = 0;
  bip->ref_counter = 1;

  bip->bs = (unsigned long long*) malloc (bip->n->ints * sizeof (unsigned long long));
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = 0ULL;

  return bip;
}

void
del_bipartition (bipartition bip)
{
  if (bip) {
    if (--bip->ref_counter) return;
    if (bip->bs) free (bip->bs);
    del_bipsize (bip->n);
    free (bip);
  }
}

void
del_bipsize (bipsize n)
{
  if (n) {
    if (--n->ref_counter) return;
    free (n);
  }
}

void
bipsize_resize (bipsize n, int nbits)
{
  int i;
  n->bits = nbits;
  n->ints = nbits/BitStringSize + 1; // might be smaller than original bs size
  n->mask = 0ULL;
  for (i=0; i < nbits%BitStringSize; i++) n->mask |= (1ULL << i); /* disregard other bits */
}

void
bipartition_initialize (bipartition bip, int position)
{
  int i, j;
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = 0ULL;
  j = position%BitStringSize;
  i = position/BitStringSize;

  bip->bs[i] = (1ULL << j);
  bip->n_ones = 1;
}

void
bipartition_zero (bipartition bip)
{
  int i;
  for (i=0; i < bip->n->ints; i++) bip->bs[i] = 0ULL;
  bip->n_ones = 0;
}

void
bipartition_set (bipartition bip, int position)
{
  bipartition_set_lowlevel (bip, position/BitStringSize, position%BitStringSize);
}

void
bipartition_set_lowlevel (bipartition bip, int i, int j)
{
  if (bip->bs[i] & (1ULL << j)) return; // bit already set
  bip->bs[i] |= (1ULL << j);
  bip->n_ones++; /* doesn't work if we reduce space later (check replace_int_in_vector() ) */
}

void
bipartition_unset (bipartition bip, int position)
{
  bipartition_unset_lowlevel (bip, position/BitStringSize, position%BitStringSize);
}

void
bipartition_unset_lowlevel (bipartition bip, int i, int j)
{
  if (!(bip->bs[i] & (1ULL << j))) return; // bit already unset
  bip->bs[i] &= ~(1ULL << j);
  bip->n_ones--;
}

void
bipartition_copy (bipartition to, const bipartition from)
{
  int i;
  for (i=0; i < to->n->ints; i++) to->bs[i] = from->bs[i];
  to->n_ones = from->n_ones;
}

void
bipartition_OR (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] | b2->bs[i];
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = b1->n_ones + b2->n_ones; // works on topologies where b1 and b2 are disjoint
}

void
bipartition_AND (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] & b2->bs[i];
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = 0;// update_count = false should be used only when you don't care about this value (temp var)
}

void
bipartition_ANDNOT (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] & (~b2->bs[i]);
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = 0;// update_count = false should be used only when you don't care about this value (temp var)
}

void
bipartition_XOR (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] ^ b2->bs[i];
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = 0;// update_count = false should be used only when you don't care about this value (temp var)
}

void
bipartition_XORNOT (bipartition result, const bipartition b1, const bipartition b2, bool update_count)
{ /* equivalent to XOR followed by NOT */
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = b1->bs[i] ^ (~b2->bs[i]);
  result->bs[i-1] &= b1->n->mask; /* do not change last bits (do not belong to bipartition) */
  if (update_count) bipartition_count_n_ones (result);
  else result->n_ones = 0;// update_count = false should be used only when you don't care about this value (temp var)
}

void
bipartition_NOT (bipartition result, const bipartition bip)
{
  int i;
  for (i=0; i < result->n->ints; i++) result->bs[i] = ~bip->bs[i];
  result->bs[i-1] &= bip->n->mask; /* do not invert last bits (do not belong to bipartition) */
  result->n_ones = bip->n->bits - bip->n_ones;
}

void
bipartition_count_n_ones (const bipartition bip)
{
  int i;
  unsigned long long j;
  bip->n_ones = 0;
/* // Naive approach
  for (i=0; i < bip->n_ints - 1; i++) for (j=0; j < BitStringSize; j++) bip->n_ones += ((bip->bs[i] >> j) & 1ULL);
  for (j=0; j < bip->n_bits%BitStringSize; j++) bip->n_ones += ((bip->bs[i] >> j) & 1ULL);
 */
  // clear the least significant bit set per iteration (Peter Wegner in CACM 3 (1960), 322, mentioned in K&R)
  for (i=0; i < bip->n->ints; i++) for (j = bip->bs[i]; j; bip->n_ones++) j &= j - 1ULL;
}

bool
bipartition_is_equal (const bipartition b1, const bipartition b2)
{
  int i;
  if (b1->n_ones  != b2->n_ones)  return false;
  if (b1->n->ints != b2->n->ints) return false;
  for (i=0; i < b1->n->ints - 1; i++) if (b1->bs[i] != b2->bs[i]) return false;
  b1->bs[i] &= b1->n->mask; b2->bs[i] &= b2->n->mask; /* apply mask before comparing last elems */
  if (b1->bs[i] != b2->bs[i]) return false;
  return true;
}

bool
bipartition_is_equal_bothsides (const bipartition b1, const bipartition b2)
{
  int i;
  bool equal = true;
  for (i=0; (i < b1->n->ints - 1) && (equal); i++) if (b1->bs[i] != b2->bs[i]) equal = false;
  if ((equal) && ((b1->bs[i] & b1->n->mask) != (b2->bs[i] & b2->n->mask))) equal = false;
  if (equal) return true; /* the biparitions are already the same, without flipping the bits */
  /* now we compare one bipartition with the complement of the other */
  for (i=0; (i < b1->n->ints - 1); i++) if (b1->bs[i] != ~b2->bs[i]) return false;
  if ((b1->bs[i] & b1->n->mask) != ((~b2->bs[i]) & b2->n->mask)) return false;
  return true; /* they are the exact complement of one another */
}

bool
bipartition_is_larger (const bipartition b1, const bipartition b2)
{
  int i;
  if (b1->n_ones > b2->n_ones) return true;
  if (b1->n_ones < b2->n_ones) return false;

  for (i = b1->n->ints - 1; (i >= 0) && (b1->bs[i] == b2->bs[i]); i--); /* find position of distinct bipartition elem*/

  if (i < 0) return false; /* identical bipartitions */
  if (b1->bs[i] > b2->bs[i]) return true;
  else return false;
}

void
bipartition_flip_to_smaller_set (bipartition bip)
{
  int i = bip->n->ints - 1; /* most significant position -- consistent with is_larger() above, using OLD algo below */
  if ((2 * bip->n_ones) < bip->n->bits) return; /* it is already the smaller set */
  /* OLD always x is different from ~x, so we just look at last element ("largest digits of number") */
  // if (((2 * bip->n_ones) == bip->n->bits) && (bip->bs[i] < (bip->n->mask & ~bip->bs[i]))) return;
  /* NEW: resolve ties by always showing the same "side" of bipartition, that is, the one having an arbitrary leaf (first one, in our case) */
  if (((2 * bip->n_ones) == bip->n->bits) && (bip->bs[0] & 1ULL)) return;

  for (i=0; i < bip->n->ints; i++) bip->bs[i] = ~bip->bs[i]; /* like bipartition_NOT() */
  bip->bs[i-1] &= bip->n->mask; /* do not invert last bits (do not belong to bipartition) */
  bip->n_ones = bip->n->bits - bip->n_ones;
  return;
}

bool
bipartition_is_bit_set (const bipartition bip, int position)
{
  if (bip->bs[(int)(position/BitStringSize)] & (1ULL << (int)(position%BitStringSize))) return true;
  return false;
}

bool
bipartition_contains_bits (const bipartition b1, const bipartition b2)
{ /* generalization of bipartition_is_bit_set(); b1 contains or not b2 */
  int i;
  if (b1->n_ones < b2->n_ones) return false;
  for (i=0; i < b1->n->ints; i++) if ((b2->bs[i]) && (b2->bs[i] != (b1->bs[i] & b2->bs[i]))) return false;
  return true;
}

void
bipartition_to_int_vector (const bipartition b, int *id, int vecsize)
{
  int i, j, k = 0;
  for (i=0; i < b->n->ints; i++) for (j=0; (j < BitStringSize) && (k < vecsize); j++) if ( ((b->bs[i] >> j) & 1ULL) ) id[k++] = i * BitStringSize + j;
}

/*
void
bipartition_print_to_stdout (const bipartition b1)
{
  int i, j;
  for (i = 0; i < b1->n->ints - 1; i++) {
    for (j = 0; j < BitStringSize; j++) printf ("%d", (int)((b1->bs[i] >> j) & 1ULL));
    printf (".");
  }
  for (j = 0; j < b1->n->bits%BitStringSize; j++) printf ("%d", (int)((b1->bs[i] >> j) & 1ULL));
  printf ("[%d] ", b1->n_ones);
}
*/

void
bipartition_replace_bit_in_vector (bipartition *bvec, int n_b, int to, int from, bool reduce)
{ /* copy info from position "from" to position "to" */
  int k, j = from%BitStringSize, i = from/BitStringSize, j2 = to%BitStringSize, i2 = to/BitStringSize;

  /* boolean "reduce" means that bitstring space will be reduced (last bits will be removed), therefore the update of
   * n_ones is different from default bipartition_set() behaviour: it's not an extra "1" (that is, one that did not
   * contribute to n_ones), but an existing "1" that change places. Schematically:
   * from -> to | normal n_ones count | when bitstring is reduced afterwards
   * 0    -> 0  |     0               |  0
   * 0    -> 1  |    -1               | -1
   * 1    -> 0  |    +1               |  0  (since it's a leaf that belonged to position "from" and now is on position "to")
   * 1    -> 1  |     0               | -1  (in fact one of the two "1"s dissapeared after reducing the bitstring
   * (the above description is outdated since I rewrote by hand the bit functions -- observe how we must erase 1 values from "from") */
  if (reduce) for (k = 0; k < n_b; k++) { // copy 0 or 1 values, erasing "from" values to avoid problems after reducing space (hanging 1s out of range)
    if      ( ((bvec[k]->bs[i] >> j) & 1ULL) && ((bvec[k]->bs[i2] >> j2) & 1ULL) )  { bvec[k]->n_ones--; bvec[k]->bs[i] &= ~(1ULL << j); }
    else if ( ((bvec[k]->bs[i] >> j) & 1ULL) && !((bvec[k]->bs[i2] >> j2) & 1ULL) ) { bvec[k]->bs[i2] |=  (1ULL << j2); bvec[k]->bs[i] &= ~(1ULL << j); }
    else if ( !((bvec[k]->bs[i] >> j) & 1ULL) && ((bvec[k]->bs[i2] >> j2) & 1ULL) ) { bvec[k]->bs[i2] &= ~(1ULL << j2); bvec[k]->n_ones--; }
    /* else do nothing (from zero to zero) */
  }

  else for (k = 0; k < n_b; k++) { // copy 0 or 1 values
    if ( ((bvec[k]->bs[i] >> j) & 1ULL) ) bipartition_set_lowlevel   (bvec[k], i2, j2); // will check if n_ones change or not
    else                                 bipartition_unset_lowlevel (bvec[k], i2, j2);
  }
}

void
bipartition_resize_vector (bipartition *bvec, int n_b)
{
  int k, i = bvec[0]->n->ints - 1;
  for (k = 0; k < n_b; k++) { bvec[k]->bs[i] &= bvec[0]->n->mask; bipartition_count_n_ones (bvec[k]); }
}


