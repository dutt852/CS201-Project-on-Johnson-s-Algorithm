#include <bits/stdc++.h>
using namespace std;

struct binomial_node { //node of binomial tree
    pair<int,int> p; //p.first=>distance , p.second=>vertex number
	int rank;
    binomial_node *child,*parent,*sibling;
};

struct fib_node {  //node of fibonacci heap
	pair<int,int> p; //p.first=>distance , p.second=>vertex number
	int rank,marked;
	fib_node *child,*parent,*right,*left;
};

binomial_node* add_new_node_binomial(int distance, int vertex) {  //adding new node in binomial heap
    binomial_node* new_node=new binomial_node;
    new_node->p.first=distance;
    new_node->p.second=vertex;
    new_node->parent=NULL;
    new_node->child=NULL;
    new_node->sibling=NULL;
    new_node->rank=0;
    return new_node;
}

vector<binomial_node*> Union(vector<binomial_node*> h1, vector<binomial_node*> h2) {  //union of two binomial heaps h1 and h2
    vector<binomial_node*> new_heap,new_heap2;
    int i1=0,i2=0,i3=0;
	//combining h1 and h2 as new_heap (non-decreasing order of ranks):-
	while(i1<h1.size() && i2<h2.size()) {
		if(((h1[i1])->rank)>((h2[i2])->rank)) {
			new_heap.push_back(h2[i2]);
			i2++;
		}
		else if(((h1[i1])->rank)<=((h2[i2])->rank)) {
			new_heap.push_back(h1[i1]);
			i1++;
		}
	}
	if(i2==h2.size()) {
		while(i1<h1.size()) {
			new_heap.push_back(h1[i1]);
			i1++;
		}
	}
	if(i1==h1.size()) {
		while(i2<h2.size()) {
			new_heap.push_back(h2[i2]);
			i2++;
		}
	}
    //re-arranging new_heap so that no two binomial trees have same rank and ranks in increasing order:-
    if(new_heap.size()<2) {
    	return new_heap;
	}
	else {
		i1=0;
		i2=1;
		i3=2;
		if(new_heap.size()==2) {
			i3--;
		}
		while(i1<new_heap.size()) {
			if(i2==new_heap.size()) {
				i1++;
			}
			else if(i3<new_heap.size() && ((new_heap[i1]->rank)==(new_heap[i2]->rank) && (new_heap[i2]->rank)==(new_heap[i3]->rank))) {
				i3++;
				i2++;
				i1++;
			}
			else if((new_heap[i2]->rank)!=(new_heap[i1]->rank)) {
				i3++;
				i2++;
				i1++;
			}
			else if((new_heap[i1]->rank)==(new_heap[i2]->rank)) {  //this case requires merging two binomial trees
					i3++;
				//to merge two binomial trees and making tree with larger distance as a child of tree with smaller distance:-
				if((new_heap[i1]->p.first)<(new_heap[i2]->p.first)) {
			    	new_heap[i2]->sibling=new_heap[i1]->child;
					new_heap[i2]->parent=new_heap[i1];
					new_heap[i1]->rank=new_heap[i1]->rank+1;
			    	new_heap[i1]->child=new_heap[i2];
			    	new_heap.erase(new_heap.begin()+i2);
				}
				else {
					new_heap[i1]->sibling=new_heap[i2]->child;
					new_heap[i1]->parent=new_heap[i2];
					new_heap[i2]->rank=new_heap[i2]->rank+1;
			    	new_heap[i2]->child=new_heap[i1];
			    	new_heap.erase(new_heap.begin()+i1);
				}
			}
		}
		return new_heap;
	}
}

vector<binomial_node*> insert_new_val_binomial(vector<binomial_node*> heap, int distance, int vertex, map<int,binomial_node*> &dict) {  //inserting a node in binomial heap
    binomial_node* new_node=add_new_node_binomial(distance,vertex);
    dict[vertex]=new_node;
    vector<binomial_node*> new_heap;
    new_heap.push_back(new_node);
    new_heap=Union(heap,new_heap);
    return new_heap;
}

int find_min(vector<binomial_node*> heap) {  //returns index corresponding to min-distance in binomial heap
	int i=0,min_index=0,min_dist=heap[0]->p.first;
	for(i=0;i<heap.size();++i) {
		if((heap[i]->p.first)<min_dist) {
			min_dist=heap[i]->p.first;
			min_index=i;
		}
	}
	return min_index;
}

vector<binomial_node*> extract_min_binomial(vector<binomial_node*> heap) {  //node corresponding to minimum value is deleted from binomial heap
	vector<binomial_node*> new_heap,new_heap2;
	binomial_node *temp,*temp2;
	int min_index,i=0;
	min_index=find_min(heap);
	temp=heap[min_index]->child;
	for(i=0;i<heap.size();++i) {
		if(i!=min_index) {
			new_heap.push_back(heap[i]);
		}
	}
	while(temp!=NULL) {
		binomial_node *next;
		next=temp->sibling;
		temp->sibling=NULL;
		new_heap2.insert(new_heap2.begin(),temp);
		temp=next;
	}
	new_heap=Union(new_heap,new_heap2);
	return new_heap;
}

void decrease_key_binomial(map<int,binomial_node*> &dict, int vertex, int new_distance) {  //decrease key in binomial heap
    binomial_node* picked_node=dict[vertex];
    picked_node->p.first=new_distance;
    while((picked_node->parent)!=NULL) {
        if((picked_node->parent->p.first)>(picked_node->p.first)) {
            int a=picked_node->p.first;
            int b=picked_node->p.second;
            picked_node->p.first=picked_node->parent->p.first;
            picked_node->p.second=picked_node->parent->p.second;
            picked_node->parent->p.first=a;
            picked_node->parent->p.second=b;
            dict[picked_node->parent->p.second]=picked_node->parent;
            dict[picked_node->p.second]=picked_node;
            picked_node=picked_node->parent;
        }
        else {
            break;
        }
    }
}

fib_node* add_new_node_fib(int distance, int vertex) {  //adding new node in fibonacci heap
	fib_node* new_node=new fib_node;
	new_node->p.first=distance;
    new_node->p.second=vertex;
    new_node->right=new_node;
    new_node->left=new_node;
    new_node->child=NULL;
    new_node->parent=NULL;
    new_node->marked=0;
    new_node->rank=0;
    return new_node;
}

void insert_new_val_fib(fib_node** head_ref, int distance, int vertex, map<int,fib_node*> &dict, int* total_nodes) {  //insertion in fibonacci heap
	fib_node* new_node=add_new_node_fib(distance,vertex);
	dict[vertex]=new_node;
	(*total_nodes)++;
	if((*head_ref)==NULL) {  //heap was empty
		*head_ref=new_node;
	}
	else if((*head_ref)!=NULL) {
		new_node->left=(*head_ref)->left;
		((*head_ref)->left)->right=new_node;
		new_node->right=(*head_ref);
		(*head_ref)->left=new_node;
		if((new_node->p.first)<((*head_ref)->p.first)) {
			*head_ref=new_node;  //head points to minimum distance fib_node
		}
	}
}

fib_node* extract_min_fib(fib_node** head_ref, int* total_nodes) {  //returns extracted node from fibonacci heap
	fib_node *x=(*head_ref),*child1;
	int flag=0;
	if(x==NULL) {  //this case will never occur actually in my implementation of Johnson algo.
		return x;
	}
	(*total_nodes)--;
	child1=x->child;
	if(child1!=NULL) {
		child1->left->right=x->right;
		x->right->left=child1->left;
		x->right=child1;
		child1->left=x;
	}
	//deleting x from root list:-
	if(x==(x->left)) {  //single node
		*head_ref=NULL;
		x->child=NULL;
		return x;
	}
	else {
		(x->right)->left=x->left;
		(x->left)->right=x->right;
		*head_ref=x->right;
		x->left=x;
		x->right=x;
		x->child=NULL;
		//consolidating root list after extracting minimum:-
		int max_degree=int(log2(*total_nodes))+2,j,degree;
		vector<fib_node*> v;  //to map ranks of nodes of root list
		fib_node *ptr1,*ptr2,*ptr3,*ptr4,*temp_ptr,*ptr5;
		fib_node *temp=*head_ref;
		ptr1=(*head_ref);
		int cntnodesnow=0;
		do{
			if((temp->p.first)>(ptr1->p.first))
				temp=ptr1;
			cntnodesnow++;
			ptr1=ptr1->right;
		}
		while(ptr1!=(*head_ref));
		for(j=0;j<max_degree+1;++j) {
			v.push_back(NULL);
		}
		int cnttemp=0;
		*head_ref=temp;
		ptr1=temp->right;
		while(cnttemp<cntnodesnow){
			cnttemp++;
			ptr2=ptr1->right;
			ptr3=ptr1;
			degree=ptr3->rank;
			while(v[degree]!=NULL) {
				ptr4=v[degree];
				if((ptr4->p.first)<(ptr3->p.first)) {
					swap(ptr3,ptr4);
				}
				//link ptr3 and ptr4 as parent-child:-
				(ptr4->right)->left=ptr4->left;
				(ptr4->left)->right=ptr4->right;
				ptr4->left=ptr4;
				ptr4->right=ptr4;
				ptr4->parent=ptr3;
				(ptr3->rank)++;
				if(ptr3->child!=NULL) {
					temp_ptr=(ptr3->child)->right;
					(ptr3->child)->right=ptr4;
					ptr4->left=ptr3->child;
					temp_ptr->left=ptr4;
					ptr4->right=temp_ptr;
				}
				else if(ptr3->child==NULL) {
					ptr4->left=ptr4;
					ptr4->right=ptr4;
					ptr3->child=ptr4;
				}  //(linking completed)
				v[degree]=NULL;
				degree++;
			}
			v[degree]=ptr3;
			v[degree]->right=v[degree];
			v[degree]->left=v[degree];
			ptr1=ptr2;
		}
		*head_ref==NULL;
		for(j=0;j<max_degree+1;++j) {
			if(v[j]!=NULL) {
				if((*head_ref)!=NULL) {
					ptr5=(*head_ref)->left;
					ptr5->right=v[j];
					v[j]->right=(*head_ref);
					v[j]->left=ptr5;
					(*head_ref)->left=v[j];
					if((v[j]->p.first)<((*head_ref)->p.first)) {  //updating head/minimum
						(*head_ref)=v[j];
					}
				}
				else if((*head_ref)==NULL) {
					(*head_ref)=v[j];
				}
			}
		}  //(consolidation completed)
		return x;
	}
}

void cutting(fib_node** head_ref, fib_node* node1, fib_node* node2) {  //node2 is actually parent of node1 in fibonacci heap
	(node2->rank)--;  //as child will be removed
	//removal of node1 from list of children:-
	if(((node2->child)->right)==(node2->child)) {  //single child case
		node2->child=NULL;
	}
	else {
		if(node2->child==node1) {
			node2->child=node1->left;
		}
		(node1->right)->left=node1->left;
		(node1->left)->right=node1->right;
	}
	//now add node1 into root list:-
	node1->marked=0;
	node1->parent=NULL;
	((*head_ref)->left)->right=node1;
	node1->right=(*head_ref);
	node1->left=((*head_ref)->left);
	(*head_ref)->left=node1;
	if(((*head_ref)->p.first)>(node1->p.first))
		*head_ref=node1;
}

void recursive_cutting(fib_node** head_ref, fib_node* node1) {
	if(node1!=NULL) {
		if((node1->parent)!=NULL) {
			if((node1->marked)==0) {  //recursion is stopped if parent is unmarked
				node1->marked=1;
			}
			else {  //recursion case
				cutting(head_ref,node1,node1->parent);
				recursive_cutting(head_ref,node1->parent);
			}
		}
		else
			node1->marked=0;
	}
}

void decrease_key_fib(fib_node** head_ref, map<int,fib_node*> &dict, int vertex, int new_distance) {  //decrease key in fibonacci heap
	dict[vertex]->p.first=new_distance;
	if(((dict[vertex]->parent)!=NULL) && (((dict[vertex]->parent)->p.first)>(dict[vertex]->p.first))) {
		cutting(head_ref,dict[vertex],dict[vertex]->parent);
		recursive_cutting(head_ref,dict[vertex]->parent);
	}
	//updating head/minimum (if required):-
	if(((*head_ref)->p.first)>(dict[vertex]->p.first)) {
		(*head_ref)=dict[vertex];
	}
}

void percolate_down(vector<pair<int,int>> &arr, int i, map<int,int> &dict) {  //in case of deletion in binary heap
	int minimum=i,left=((2*i)+1),right=((2*i)+2),temp_v,temp_wt;  //indices
	if(left<arr.size()) {  //left child exists
		if(arr[left].first<arr[minimum].first) {
			minimum=left;
		}
	}
	if(right<arr.size()) {  //right child exists
		if(arr[right].first<arr[minimum].first) {
			minimum=right;
		}
	}
	if(minimum!=i) {
		//swap arr[minimum] and arr[i] and update dict accordingly:
		dict[arr[i].second]=minimum;
		dict[arr[minimum].second]=i;
		temp_v=arr[i].second;
		temp_wt=arr[i].first;
		arr[i].first=arr[minimum].first;
		arr[i].second=arr[minimum].second;
		arr[minimum].first=temp_wt;
		arr[minimum].second=temp_v;
		//and repeat same procedure for sub-tree using recursion
		percolate_down(arr,minimum,dict);
	}
}

void delete_min_val_binary(vector<pair<int,int>> &arr, map<int,int> &dict) {  //at root of binary heap
	if(arr.size()>0) {
		int temp_v,temp_wt;
		//swap root with last node and update dict accordingly:
		dict[arr[0].second]=arr.size()-1;
		dict[arr[arr.size()-1].second]=0;
		temp_v=arr[0].second;
		temp_wt=arr[0].first;
		arr[0].first=arr[arr.size()-1].first;
		arr[0].second=arr[arr.size()-1].second;
		arr[arr.size()-1].first=temp_wt;
		arr[arr.size()-1].second=temp_v;
		arr.pop_back();  //reducing heap size
		percolate_down(arr,0,dict);
	}
}

void percolate_up(vector<pair<int,int>> &arr, int i, map<int,int> &dict) {  //in case of insertion or decrease key in binary heap
	int parent,temp_v,temp_wt;  //index
	if(i>0) {
		parent=(int)((i-1)/2);
		if(arr[parent].first>arr[i].first) {
			//swap and update dict accordingly:-
			dict[arr[i].second]=parent;
			dict[arr[parent].second]=i;
			temp_v=arr[i].second;
			temp_wt=arr[i].first;
			arr[i].first=arr[parent].first;
			arr[i].second=arr[parent].second;
			arr[parent].first=temp_wt;
			arr[parent].second=temp_v;
			//and repeat same procedure for sub-tree using recursion:-
			percolate_up(arr,parent,dict);
		}
	}
}

void insert_new_val_binary(vector<pair<int,int>> &arr, int new_wt, int new_v, map<int,int> &dict) {  //insertion in binary heap
	arr.push_back(make_pair(new_wt,new_v));
	dict[new_v]=arr.size()-1;
	percolate_up(arr,arr.size()-1,dict);
}

void decrease_key_binary(vector<pair<int,int>> &arr, int index, int new_wt, map<int,int> &dict) {  //decrease key in binary heap
	arr[index].first=new_wt;
	percolate_up(arr,index,dict);
}

void bellman_ford(vector<pair<int,int>> adj[], int n, int h[], int* is_applicable) {
	int i,j,k,v,wt,change,s=n;  //s->source
	for(i=0;i<n+1;++i) {
		h[i]=999999;
	}
	h[s]=0;
	for(k=0;k<n+1;++k) {  //n+1 times relax for all edges ((n+1)th time is to check if there is a negative cycle)
		change=0;
		for(i=0;i<n+1;++i) {
			for(j=0;j<adj[i].size();++j) {
				v=adj[i][j].first;
				wt=adj[i][j].second;
				if(h[v]>(h[i]+wt)) {  //relax operation
					h[v]=h[i]+wt;
					change=1;
				}
			}
		}
		if(change==0) {
			break;
		}
		if(k==n && change!=0) {
			*is_applicable=0;
		}
	}
}

void dijkstra_binomial(vector<pair<int,int>> adj[], int n, int h[]) {  //binomial-heap based implementation
	int i,j,dist[n],u,v,wt,s;
	for(j=0;j<n;++j) {
		int visited[n]={n*0};
		s=j;  //source
		vector<binomial_node*> heap;  //binomial heap
		map<int,binomial_node*> dict;  //for mapping vertices to corresponding heap nodes
		for(i=0;i<n;++i) {
			dist[i]=999999;
		}
		for(int i=0;i<n;i++){
    		dict[i]=NULL;
		}
		dist[s]=0;
		heap=insert_new_val_binomial(heap,dist[s],s,dict);
		visited[s]=1;
		while(heap.size()>0) {
			u=heap[find_min(heap)]->p.second;
			heap=extract_min_binomial(heap);
			for(i=0;i<adj[u].size();++i) {
				v=adj[u][i].first;
				wt=adj[u][i].second;
				if(dist[v]>(dist[u]+wt)) {  //relax operation
					dist[v]=dist[u]+wt;
					if(visited[v]==0) {
						heap=insert_new_val_binomial(heap,dist[v],v,dict);
						visited[v]=1;
					}
					else {
						decrease_key_binomial(dict,v,dist[v]);
					}
				}
			}
		}
		for(i=0;i<n;++i) {
			if(dist[i]==999999) {
				cout <<dist[i]<<" ";
			}
			else {
				cout <<dist[i]+h[i]-h[s]<<" ";
			}
		}
		cout <<"\n";
	}
}

void dijkstra_binary(vector<pair<int,int>> adj[], int n, int h[]) {  //binary-heap based implementation
	int i,j,dist[n],u,v,wt,s;
	for(j=0;j<n;++j) {
		int visited[n]={n*0};
		s=j;  //source
		vector<pair<int,int>> arr;  //binary heap
		map<int,int> dict;  //for mapping vertices to corresponding heap indices
		for(i=0;i<n;++i) {
			dist[i]=999999;
		}
		dist[s]=0;
		insert_new_val_binary(arr,dist[s],s,dict);
		visited[s]=1;
		while(arr.size()>0) {
			u=arr[0].second;
			delete_min_val_binary(arr,dict);
			for(i=0;i<adj[u].size();++i) {
				v=adj[u][i].first;
				wt=adj[u][i].second;
				if(dist[v]>(dist[u]+wt)) {  //relax operation
					dist[v]=dist[u]+wt;
					if(visited[v]==0) {
						insert_new_val_binary(arr,dist[v],v,dict);
						visited[v]=1;
					}
					else {
						decrease_key_binary(arr,dict[v],dist[v],dict);
					}
				}
			}
		}
		for(i=0;i<n;++i) {
			if(dist[i]==999999) {
				cout <<dist[i]<<" ";
			}
			else {
				cout <<dist[i]+h[i]-h[s]<<" ";
			}
		}
		cout <<"\n";
	}
}

void dijkstra_array(vector<pair<int,int>> adj[], int n, int h[]) {  //array based implementation
	int i,j,k,dist[n],u,v,wt,s,min_val,min_index,temp_v,temp_wt;
	for(j=0;j<n;++j) {
		int visited[n]={n*0};
		s=j;  //source
		vector<pair<int,int>> arr;
		for(i=0;i<n;++i) {
			dist[i]=999999;
		}
		dist[s]=0;
		arr.push_back(make_pair(dist[s],s));
		visited[s]=1;
		while(arr.size()>0) {
			if(arr.size()>1) {  //finding minimum in array
				min_val=arr[0].first;
				min_index=0;
				for(i=0;i<arr.size();++i) {
					if(arr[i].first<min_val) {
						min_val=arr[i].first;
						min_index=i;
					}
				}
				temp_v=arr[min_index].second;
				temp_wt=arr[min_index].first;
				arr[min_index].first=arr[arr.size()-1].first;
				arr[min_index].second=arr[arr.size()-1].second;
				arr[arr.size()-1].first=temp_wt;
				arr[arr.size()-1].second=temp_v;
			}
			u=arr[arr.size()-1].second;
			arr.pop_back();
			for(i=0;i<adj[u].size();++i) {
				v=adj[u][i].first;
				wt=adj[u][i].second;
				if(dist[v]>(dist[u]+wt)) {  //relax operation
					dist[v]=dist[u]+wt;
					if(visited[v]==0) {
						arr.push_back(make_pair(dist[v],v));
						visited[v]=1;
					}
					else {
						for(k=0;k<arr.size();++k) {
							if(arr[k].second==v) {  //decrease key
								arr[k].first=dist[v];
								break;
							}
						}
					}
				}
			}
		}
		for(i=0;i<n;++i) {
			if(dist[i]==999999) {
				cout <<dist[i]<<" ";
			}
			else {
				cout <<dist[i]+h[i]-h[s]<<" ";
			}
		}
		cout <<"\n";
	}
}

void dijkstra_fib(vector<pair<int,int>> adj[], int n, int h[]) {  //fibonacci-heap based implementation
	int i,j,k,dist[n],u,v,wt,s;
	for(j=0;j<n;++j) {
		int visited[n]={n*0};
		s=j;  //source
		fib_node* head=NULL;  //fibonacci heap's head pointer or minimum pointer
		fib_node* extracted_node;
		int total_nodes=0;  //to count total nodes in fibonacci heap
		map<int,fib_node*> dict;  //for mapping vertices to corresponding heap nodes
		for(i=0;i<n;++i) {
			dist[i]=999999;
		}
		for(int i=0;i<n;i++){
    		dict[i]=NULL;
		}
		dist[s]=0;
		insert_new_val_fib(&head,dist[s],s,dict,&total_nodes);
		visited[s]=1;
		while(head!=NULL) {
			extracted_node=extract_min_fib(&head,&total_nodes);
			u=extracted_node->p.second;
			for(i=0;i<adj[u].size();++i) {
				v=adj[u][i].first;
				wt=adj[u][i].second;
				if(dist[v]>(dist[u]+wt)) {  //relax operation
					dist[v]=dist[u]+wt;
					if(visited[v]==0) {
						insert_new_val_fib(&head,dist[v],v,dict,&total_nodes);
						visited[v]=1;
					}
					else {
						decrease_key_fib(&head,dict,v,dist[v]);
					}
				}
			}
		}
		for(i=0;i<n;++i) {
			if(dist[i]==999999) {
				cout <<dist[i]<<" ";
			}
			else {
				cout <<dist[i]+h[i]-h[s]<<" ";
			}
		}
		cout <<"\n";
	}
}

int main(int argc, char* argv[]) {
	int T;
	cin>>T;
	vector <double> time_arr;
	while(T--) {
		clock_t start,end;
		start=clock();
		int n,d;
		cin >>n>>d;
		int i,j,val,h[n+1]={(n+1)*0},neg_edge=0,is_applicable=1;
		vector<pair<int,int>> adj[n+1];
		for(i=0;i<n;++i) {
			for(j=0;j<n;++j) {
				cin >>val;  //edge-weight
				if(val<0) {
					neg_edge=1;
				}
				if(val!=999999 && i!=j) {
					adj[i].push_back(make_pair(j,val));
				}
			}
		}
		for(j=0;j<n;++j) {
			adj[n].push_back(make_pair(j,0));  //new vertex with zero-weight outgoing edges to all other vertices in original graph
		}
		if(neg_edge==1 && d==0) {
			is_applicable=0;
		}
		if(neg_edge==1 && d!=0) {  //if graph is directed with atleast one negative edge
			bellman_ford(adj,n,h,&is_applicable);
			if(is_applicable!=0) {  //update edge-weights
				for(i=0;i<n;++i) {
					for(j=0;j<adj[i].size();++j) {
						adj[i][j].second=adj[i][j].second+h[i]-h[adj[i][j].first];
					}
				}
			}
		}
		if(is_applicable==0) {
			cout <<-1<<"\n";
		}
		else {
			//now we have to apply n times dijkstra on updated graph
			if(argc==1) {
				dijkstra_fib(adj,n,h);
			}
			else {
				int ch;
				ch=atoi(argv[1]);
				if(ch==1) {
					dijkstra_array(adj,n,h);
				}
				else if(ch==2) {
					dijkstra_binary(adj,n,h);
				}
				else if(ch==3) {
					dijkstra_binomial(adj,n,h);
				}
				else {
					dijkstra_fib(adj,n,h);
				}
			}
		}
		end=clock();
		//total time taken by one test case:-
		double TakenTime=((double)(end-start))/CLOCKS_PER_SEC;
		time_arr.push_back(TakenTime);
	}
	for(int I=0;I<time_arr.size();++I) {
		cout <<time_arr[I]<<" ";
	}
	cout <<"\n";
	return 0;
}
