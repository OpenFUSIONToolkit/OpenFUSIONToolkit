/*-----------------------------------------------------------------------------
* Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
*
* SPDX-License-Identifier: LGPL-3.0-only
*------------------------------------------------------------------------------
* Sorting and search functions for the Open FUSION Toolkit
*----------------------------------------------------------------------------*/
#include <algorithm>
#include <vector>

template <class T1, class T2>
void int_sort(T1* ia,T2 n) {
	std::sort(ia,ia+n);
}

template <class T1, class T2>
T1 int_search(T1* list,T2 length, T1 item) {
	T2 ll;
	T2 lu;
	T2 lm;
	T2 result;
	ll=0;
	lu=length-1;
	if(length==0) {
		result=0;
		return result;
	}
	while(lu-ll>1) {
		lm=(ll+lu)/2;
		if(list[lm]-item<0)
			ll=lm;
		else if (list[lm]-item==0) {
			result=lm+1;
			return result;
		} else if (list[lm]-item>0)
			lu=lm;
	}
	if(list[ll]==item){
		result=ll+1;
		return result;
	} else if (list[lu]==item) {
		result=lu+1;
		return result;
	} else {
		result=0;
		return result;
	}
}

template <typename Container>
struct compare_indirect_index
{
	const Container& container;
	compare_indirect_index( const Container& container ): container( container ) { }
	bool operator () ( size_t lindex, size_t rindex ) const
	{
		return container[ lindex ] < container[ rindex ];
	}
};

template <class T1, class T2>
void sort_array(T1* col, T2* ind, T2 n) {
	std::vector<size_t> indices(n);
	for(int i = 0; i<n;i++)
	{
		  indices[i] = i;
	}
	std::sort( indices.begin(), indices.end(), compare_indirect_index <T1*> ( col ));
	std::vector<T1> col_copy(n);
	std::vector<T2> ind_copy(n);
	for(int i = 0; i<n;i++)
	{
		col_copy[i] = col[i];
		ind_copy[i] = ind[i];
	}
	for(int i = 0; i<n;i++)
	{
		col[i] = col_copy[indices[i]];
		ind[i] = ind_copy[indices[i]];
	}
}

template <class T1, class T2>
void sort_matrix(T1* mat,T2* ind, T2 n) {
	T1 store_col1;
	T1 store_col2;
	T2 k;
	// Order row entries (min,max)
	for (T2 j=0; j<n; j++) {
		if(mat[j*2]>mat[j*2+1]) { // Row order is (max,min)
			store_col1=mat[j*2]; store_col2=mat[j*2+1];
			mat[j*2]=store_col2; mat[j*2+1]=store_col1;
			ind[j]=-ind[j];
		}
	}
	for (T2 j=0; j<n; j++) {
		for (T2 i=n-1; i>j; i--) {
			if(mat[(i-1)*2]>mat[i*2] || (mat[(i-1)*2]==mat[i*2] && mat[(i-1)*2+1]>mat[i*2+1])) { //Swap rows and indices
				store_col1=mat[i*2]; store_col2=mat[i*2+1]; k=ind[i]; // [1]
				mat[i*2]=mat[(i-1)*2]; mat[i*2+1]=mat[(i-1)*2+1]; ind[i]=ind[i-1]; // [2]
				mat[(i-1)*2]=store_col1; mat[(i-1)*2+1]=store_col2; ind[i-1]=k; // [3]
			}
		}
	}
}

extern "C" {
	void int_sort88(long* ia,long* n) { int_sort<long,long>(ia,*n);}
	void int_sort48(int* ia,long* n) { int_sort<int,long>(ia,*n);}
	void int_sort84(long* ia,int* n) { int_sort<long,int>(ia,*n);}
	void int_sort44(int* ia,int* n) { int_sort<int,int>(ia,*n);}

	void int_search88(long* list,long* length, long* item, long* result) { *result=int_search<long,long>(list,*length,*item);}
	void int_search48(int* list,long* length, int* item, long* result) { *result=int_search<int,long>(list,*length,*item);}
	void int_search84(long* list,int* length, long* item, int* result) { *result=int_search<long,int>(list,*length,*item);}
	void int_search44(int* list,int* length, int* item, int* result) { *result=int_search<int,int>(list,*length,*item);}

	void int_sort_array88(long* col,long* ind,long* n) { sort_array<long,long>(col,ind,*n);}
	void int_sort_array48(int* col,long* ind,long* n) { sort_array<int,long>(col,ind,*n);}
	void int_sort_array84(long* col,int* ind,int* n) { sort_array<long,int>(col,ind,*n);}
	void int_sort_array44(int* col,int* ind,int* n) { sort_array<int,int>(col,ind,*n);}

	void real_sort_array8(double* col,long* ind,long* n) { sort_array<double,long>(col,ind,*n);}
	void real_sort_array4(double* col,int* ind,int* n) { sort_array<double,int>(col,ind,*n);}

	void int_sort_matrix88(long* mat,long* ind,long* n) { sort_matrix<long,long>(mat,ind,*n);}
	void int_sort_matrix48(int* mat,long* ind,long* n) { sort_matrix<int,long>(mat,ind,*n);}
	void int_sort_matrix84(long* mat,int* ind,int* n) { sort_matrix<long,int>(mat,ind,*n);}
	void int_sort_matrix44(int* mat,int* ind,int* n) { sort_matrix<int,int>(mat,ind,*n);}

	void real_sort_matrix8(double* mat,long* ind,long* n) { sort_matrix<double,long>(mat,ind,*n);}
	void real_sort_matrix4(double* mat,int* ind,int* n) { sort_matrix<double,int>(mat,ind,*n);}
}
