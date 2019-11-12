#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <climits>
#include <cfloat>
#include <cassert>
#include <memory.h>
#include <sys/times.h>
#include <unistd.h>
#include <cmath>
#include <set>
#include <algorithm>
#include <utility>
#include <queue>

using namespace std;

struct Edge{
    int v1;
    int v2;
};

Edge *edge;////////////////////////////////////////

int edge_num;

//存储score为0的节点
int *t;/////////////////////////////////
int t_length;
int *t_index;/////////////////////////////////////////

typedef struct Vertex_information{
    int num_in_c;
    int config;
    char is_in_c;
    int cost;
    int score;						//加入该顶点后，自己及邻居从未支配到支配的顶点数量
    int time_stamp;
    int locked;
    int is_in_search;       //是否在搜索集
}Vertex_information;

Vertex_information *cs;

class NeighborSet
{
public:
    explicit NeighborSet(int _v)
    {
        v = _v;
        neighbor_cnt = cs[v].score;
        neighbors = new int[neighbor_cnt];
    }
    ~NeighborSet()
    {
        delete neighbors;
        neighbors = nullptr;
    }
    bool is_in_set(int _v)
    {
        for (int i = 0; i < neighbor_cnt; ++i) {
            if (neighbors[i] == _v)
                return true;
        }
        return false;
    }
    int v;
    int* neighbors;
    int neighbor_cnt;
};

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

void init();
inline int compare(int s1, int c1, int s2, int c2);
void init_best();
void add(int c, int locked_add, int init_add);
void remove(int c, int init_remove);
int in_tabu(int i);
int find_best_in_c(int allowTabu);
void uncov_r_weight_inc();
void localsearch(int maxStep);
void update_best_sol();
int check();
void super_set_reduce();
void print_reduce_graph();
/**********api start*********/
//lock vertex
void lock_vertex(int c, int locked_add);
void lock_vertex1(int v);
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
    cout << "processing " << filename << endl;
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

    for(size_t i=0;i<vertex_num;i++){
        cs[i].cost=1;
        cs[i].config=2;
        cs[i].time_stamp=1;
        cs[i].is_in_c=0;
        cs[i].num_in_c=0;
        cs[i].locked = 0;
        cs[i].is_in_search = 0;
        vertex_neightbourNum[i]=0;
        best_sol[i]=0;
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