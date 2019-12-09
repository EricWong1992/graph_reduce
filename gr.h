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
#include "ilcplex/ilocplex.h"

using namespace std;

struct Edge{
    int v1;
    int v2;
};

//顶点状态
enum State
{
    Candidate = 0,          //剩余顶点
    Fixed,                            //固定，在最优解中
    Delete,                         //可删除节点
    Forbid,                         //禁止加入最优解，但又不可删除
};

Edge *edge;////////////////////////////////////////

int edge_num;

//存储score为0的节点
int *t;/////////////////////////////////
int t_length;
int *t_index;/////////////////////////////////////////

//顶点信息
typedef struct Vertex_information{
    State state;                //顶点状态
    int num_in_c;           //被支配的次数
    int cost;                       //顶点花费
    int score;						//加入该顶点后，自己及邻居从未支配到支配的顶点数量
    int is_in_search;       //是否在搜索集
    void dump()
    {
        cout << "state:" << state << endl;
    }
}Vertex_information;

//所有顶点信息
Vertex_information *cs;

//顶点的闭邻居集合
class NeighborSet
{
public:
    explicit NeighborSet(int _v)
    {
        v = _v;
    }
    explicit NeighborSet(int _v, int n_cnt)
    {
        v = _v;
        neighbors = vector<int>(n_cnt);
    }
    ~NeighborSet()= default;
    //判断_v是否在v的闭邻居集合
    bool is_in_set(int _v)
    {
        for (int i = 0; i < getNeighborCnt(); ++i) {
            if (neighbors[i] == _v)
                return true;
        }
        return false;
    }
    //添加邻居顶点, 用于构建闭邻居集合
    void addNeighbor(int _v)
    {
        neighbors.push_back(_v);
        if (cs[_v].num_in_c == 0)
        {
            m_iValidNeighborNum++;
        }
    }
    //获取有效邻居(未支配)总数
    int getValidNeighborCnt()
    {
        return m_iValidNeighborNum;
    }
    //输出邻居总数
    int getNeighborCnt()
    {
        return neighbors.size();
    }
    //输出v及v的闭邻居集合
    void dump()
    {
        cout << v;
        for (size_t i = 0; i < neighbors.size(); i++)
        {
            cout << " " << neighbors[i];
        }
        cout << endl;
    }
    //顶点
    int v;
    //v的邻居集合
    vector<int> neighbors;
    //有效邻居总数
    int m_iValidNeighborNum = 0;
};

tms start, finish;
int start_time;
double real_time;

//程序时间限制，为0时不限制时间
double time_limit;

int *best_array;
long best_value;

int vertex_num;//顶点个数

int total_lock = 0; //fix顶点个数

//未支配顶点信息
int uncover_num;
int  *uncover_vertex;
int *uncover_vertex_index;
int uncover_length;

//剩余顶点信息
int *remain_vertex;
int remain_num;

//未使用
int *cs_vertex;
int *cs_vertex_index;
int cs_size;

//顶点二维数组，第一维顶点，第二维顶点的邻居
int** vertex;													//vertex 2d-array
int *vertex_neightbourNum;					//neighbornum
int *vertex_neightbourNum_bak;			//original neighbornum backup

int *best_sol;					//最优解，未使用
int *vertex_weight;		//顶点权值
string filename;        //实例文件名
vector<int> candidate;

void init();
//检查解是否完全支配图
int check();
//初始缩减，处理图边缘顶点
void init_reduce();
//超集缩减，处理图内部顶点
void superset_reduce();
//子集缩减
void subset_reduce();
//固定顶点
void lock_vertex(int c);
//输出信息
void print_reduce_graph_info();
void print_density(int i);
void print_degree();
void generate_reduce_graph(int step);
void add_select_vertex(vector<int>, IloNumArray, int);
//释放内存
void free_all(){
    free(vertex_neightbourNum_bak);
    free(vertex_neightbourNum);
    free(best_sol);
    free(vertex_weight);
    free(t);
    free(t_index);
    free(cs);
    free(best_array);
    free(uncover_vertex);
    free(uncover_vertex_index);
    free(remain_vertex);
    free(cs_vertex_index);
    free(cs_vertex);
    for(int i = 0; i < vertex_num; i++)
        free(vertex[i]);
}

//处理实例
int build_instance_massive(char *file_name)
{
    char line[1024];
    char tempstr1[10];
    char tempstr2[10];
    int		v1,v2;

    ifstream infile(file_name);
    filename = file_name;

    cout << "processing " << file_name << endl;
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
        cs[i].state = State::Candidate;
        cs[i].cost=1;
        cs[i].num_in_c=0;
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

//是否全部支配
bool is_all_dominated()
{
    // cout << uncover_num << endl;
    return uncover_num == 0;
}