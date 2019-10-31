#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <memory.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include <set>
#include <algorithm>
#include <utility>

using namespace std;

struct Edge{
	int v1;
	int v2;
};

Edge *edge;////////////////////////////////////////

int edge_num;

int *t;/////////////////////////////////
int t_length;
int *t_index;/////////////////////////////////////////

typedef struct Vertex_information{
	int score;
	int time_stamp;
	int config;
	char is_in_c;
	int cost;
	int num_in_c;
	int locked;
}Vertex_information;

Vertex_information *cs;

tms start, finish;
int start_time;
double real_time;

int seed;
double time_limit;
int step;
int total_step;

long BEST;

int *best_array;
long best_value;

set<int> tabu_list;

int vertex_num;//顶点个数

int *reduce;
int reduce_num;

int uncover_num;
int  *uncover_vertex;
int *uncover_vertex_index;
int uncover_length;

int *remain_vertex;
int remain_num;

int *cs_vertex;
int *cs_vertex_index;
int cs_size;

int locked_num;


//int*	vertex[MAXV];//顶点集合
int** vertex;
int *vertex_neightbourNum;
int *vertex_neightbourNum1;

int *best_sol;//最优解
int *vertex_weight;

void init();
inline int compare(int s1, int c1, int s2, int c2);
void init_best();
void add(int c, int locked_add, int init_add);
void remove(int c, int init_remove);
int in_tabu(int i);
int find_best_in_c(int allowTabu);
void uncov_r_weight_inc();
void localsearch(int maxStep);
long cost_C();
void update_best_sol();
int check();

void free_all(){

    free(vertex_neightbourNum1);
	free(vertex_neightbourNum);
	free(best_sol);
	free(vertex_weight);
	free(t);
	free(t_index);
	free(cs);
	free(best_array);
	free(reduce);
	free(uncover_vertex);
	free(uncover_vertex_index);
	free(remain_vertex);
	free(cs_vertex_index);
	free(cs_vertex);
	for(int i = 0; i < vertex_num; i++)
        free(vertex[i]);
}

int build_instance_massive(char *filename)
{
	char line[1024];
	char tempstr1[10];
	char tempstr2[10];
	int  v,e;
	char	tmp;
	int		v1,v2;

	ifstream infile(filename);
    //    if(infile==NULL) return 0;
	/*** build problem data structures of the instance ***/
	infile.getline(line,1024);
	sscanf(line, "%s %s %d %d", tempstr1,tempstr2, &vertex_num, &edge_num);

    int i,j,temp;
    edge = (Edge *)malloc(edge_num*sizeof(Edge));
    t = (int *)malloc(vertex_num*sizeof(int));

    cs = (Vertex_information *)malloc(vertex_num*sizeof(Vertex_information));
    vertex_neightbourNum = (int *)malloc(vertex_num*sizeof(int));

    vertex_neightbourNum1 = (int *)malloc(vertex_num*sizeof(int));
    vertex = (int **)malloc((edge_num*2+1000)*sizeof(int));

    best_sol=(int *)malloc(vertex_num*sizeof(int));//最优解
    vertex_weight=(int *)malloc(vertex_num*sizeof(int));//动态权重


	for(j=0;j<vertex_num;j++){
		//infile>>tempstr1>>v1>>v2;
		//temp=(j+1)%200 + 1;
		//cs[j].cost=temp;
		//v1--;
		cs[j].cost=1;
		//cs[v1].cost=v2;
		cs[j].config=2;
		cs[j].time_stamp=1;
		cs[j].is_in_c=0;
		cs[j].num_in_c=0;
		cs[j].locked = 0;
		vertex_neightbourNum[j]=0;
		best_sol[j]=0;
		t[j] = 0;
	}

	for (e=0; e<edge_num; e++)	
	{
		infile>>tempstr1>>v1>>v2;
		v1--;
		v2--;
		vertex_neightbourNum[v1]++;
		vertex_neightbourNum[v2]++;
		edge[e].v1 = v1;
		edge[e].v2 = v2;
	}
	infile.close();

	/* build v_adj and v_edges arrays */
	for (v=0; v<vertex_num; v++)
        vertex[v]=new int[vertex_neightbourNum[v]];
	for (e=0; e<edge_num; e++)
	{
		v1=edge[e].v1;
		v2=edge[e].v2;

		vertex[v1][t[v1]] = v2;
		vertex[v2][t[v2]] = v1;

		t[v1]++;
		t[v2]++;
	}

    free(edge);

    t_index = (int *)malloc(vertex_num*sizeof(int));
    reduce = (int *)malloc(vertex_num*sizeof(int));
    uncover_vertex =  (int *)malloc(vertex_num*sizeof(int));
    uncover_vertex_index=(int *)malloc(vertex_num*sizeof(int));
    remain_vertex = (int *)malloc(vertex_num*sizeof(int));
    cs_vertex = (int *)malloc(vertex_num*sizeof(int));
    cs_vertex_index = (int *)malloc(vertex_num*sizeof(int));
    best_array = (int *)malloc(vertex_num*sizeof(int));

	for(i=0;i<vertex_num;i++) 
	{
		vertex_weight[i]=1;
		//cs[i].score=vertex_neightbourNum[i]+1;//自身也要算+1
		cs[i].score = 0;
		vertex_neightbourNum1[i] = vertex_neightbourNum[i];
        }
    remain_num = vertex_num;
    uncover_num = vertex_num;

	return 1;
}



