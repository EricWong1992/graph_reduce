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
	int num_in_c;
	int config;
	char is_in_c;
	int cost;
	//以下有用
	int score;						//加入该顶点后，自己及邻居从未支配到支配的顶点数量
	int time_stamp;
	int locked;
	bool dominated;
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
int** vertex;													//vertex 2d-array
int *vertex_neightbourNum;					//neighbornum
int *vertex_neightbourNum_bak;			//original neighbornum backup

int *best_sol;					//最优解
int *vertex_weight;		//顶点权值
int *dominate_set;		//支配集

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
void graph_reduce();
void print_reduce_graph();
/**********api start*********/
//get the degree of a vertex
int getDegree(int v);
//get the vertex is locked(is in mustin)
bool getIsLocked(int v);
//get the vertex is dominated(the vertex is covered)
bool getIsDominated(int v);
//get the target_v is can dominated if add_v is added
bool isCanDominated(int add_v, int target_v);
//neighbors got dominated
void neighborGotDomimated(int v);
//get the score of adding the vertex
int getScore(int v);
//add vertex to mustin
void addVertex(int v);
/**********api end*********/

void free_all(){
    free(vertex_neightbourNum_bak);
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
	char	tmp;
	int		v1,v2;

	ifstream infile(filename);
	/*** build problem data structures of the instance ***/
	infile.getline(line,1024);
	sscanf(line, "%s %s %d %d", tempstr1,tempstr2, &vertex_num, &edge_num);

    edge = (Edge *)malloc(edge_num*sizeof(Edge));
    t = (int *)malloc(vertex_num*sizeof(int));

    cs = (Vertex_information *)malloc(vertex_num*sizeof(Vertex_information));
    vertex_neightbourNum = (int *)malloc(vertex_num*sizeof(int));

    vertex_neightbourNum_bak = (int *)malloc(vertex_num*sizeof(int));
    vertex = (int **)malloc((edge_num*2+1000)*sizeof(int));

    best_sol=(int *)malloc(vertex_num*sizeof(int));//最优解
    vertex_weight=(int *)malloc(vertex_num*sizeof(int));//动态权重
	dominate_set = (int*)malloc(vertex_num * sizeof(int));//支配集

	for(size_t i=0;i<vertex_num;i++){
		cs[i].cost=1;
		cs[i].config=2;
		cs[i].time_stamp=1;
		cs[i].is_in_c=0;
		cs[i].num_in_c=0;
		cs[i].locked = 0;
		cs[i].dominated = false;
		vertex_neightbourNum[i]=0;
		best_sol[i]=0;
		dominate_set[i] = 0;
		t[i] = 0;
	}

	for (size_t e=0; e<edge_num; e++)	
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
	for (size_t v=0; v<vertex_num; v++)
                vertex[v]=new int[vertex_neightbourNum[v]];
	for (size_t e=0; e<edge_num; e++)
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

	for(size_t i=0;i<vertex_num;i++) 
	{
		vertex_weight[i]=1;
		cs[i].score = vertex_neightbourNum[i] + 1;//自身也要算+1
		vertex_neightbourNum_bak[i] = vertex_neightbourNum[i];
        }
    remain_num = vertex_num;
    uncover_num = vertex_num;

	return 1;
}

int getDegree(int v)
{
	return vertex_neightbourNum[v];
}

bool getIsLocked(int v)
{
	return cs[v].locked == 1;
}

bool getIsDominated(int v)
{
	if (getIsLocked(v)) 
		return true;
	return cs[v].dominated;
}

bool isCanDominated(int add_v, int target_v)
{
	if (getIsDominated(target_v)) return true;
	int degree = getDegree(target_v);
	for (size_t i = 0; i < degree; i++)
	{
		if (vertex[target_v][i] == add_v)
		{
			return true;
		}
	}
	return false;
}

void neighborGotDomimated(int v)
{
	int neighbor_count = vertex_neightbourNum[v];
	for (size_t i = 0; i < neighbor_count; i++)
	{
		int v_neighbor = vertex[v][i];
		cs[v_neighbor].dominated = true;
	}
}

int getScore(int v)
{
	return cs[v].score;
}

void addVertex(int v)
{
	cs[v].locked = 1;
	dominate_set[v] = 1;
	cs[v].score = 0;
	int v_neighbor_count = vertex_neightbourNum[i];
	//修改v的邻居们的被支配状态
	for (size_t i = 0; i < v_neighbor_count; i++)
	{
		int v_neighbor = vertex[v][i];
		cs[v_neighbor].dominated = 1;
	}
	//修改v的邻居们的score值
	for (size_t i = 0; i < v_neighbor_count; i++)
	{
		int new_v = vertex[v][i];
		int new_v_neighbor_count = vertex_neightbourNum[new_v];
		int score = 0;
		for (size_t j = 0; j < new_v_neighbor_count; j++)
		{
			int new_v_neighbor = vertex[new_v][j];
			if (!getIsLocked(new_v_neighbor) && !getIsDominated(new_v_neighbor))
				score++;
		}
		cs[new_v].score = score;
	}
}