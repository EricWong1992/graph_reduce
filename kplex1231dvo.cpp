#include <iostream>
#include <cstring>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <queue>
#include <algorithm>
#include<unordered_set>
#include<sys/time.h>
#include<sys/types.h>
#include<sys/resource.h>
#include<limits.h>
// #include "Array.h"
#include "mds.h"

//  #define DEBUG

using namespace std;

#define mem_2G 536870912
#define insert_v(end, value) *(end++) = value;
#define delete_i(index, end) *index = *(--end);
#define fix(i) vertex[i].state = State::Fixed;
#define reduce(i) vertex[i].state = State::Delete;

int v_n;//顶点个数
int e_n;//边个数

int run_time = 1000;//time_limit，运行时间阈值
int remain_num;//剩余顶点数量
int uncover_num;//未支配顶点数量
int fix_num;//固定顶点数量

int is_reduce = 1;
//----------辅助的顶点邻域信息，在BB过程中记录当次迭代的G'
int *neighbor_len;//
int **neighbor;//
//int *state;//the sate of vertex
//---------------

int *temp_array;//辅助数组，用来记录程序中一些辅助信息，具体辅助信息根据需要改变
int *temp_mark;//辅助数组，用来记录程序中一些辅助信息
int *temp_index;//辅助数组，用来记录程序中一些辅助信息
int *neighbor_in_solution;//用来记录顶点与当前解中的邻居个数
//int *kv;
//int *v_weight;
Edge *edge;//用来记录边，在初始化以后因为不用到了，就释放掉
Vertex_sort *vertex_sort;//用来顶点排序
Vertex *vertex;//顶点数组

Array *crusol;//当前解数组
Array *remaining_vertex;//初始化约减后存储剩余的顶点
int para_k = 1;//kplex的参数k

Array *U;//未决策顶点集
Array *S;//当前解集，作用相当于crusol，S = crusol

int *best_solution;//用来记录最优解
int best_solution_size = 0;//最优解规模

int max_degree = -1;//最大顶点度数

double read_time, init_time, search_time;//读入时间，初始化约减时间，搜索时间
int UB = INT_MAX;//上界
int LB = 0;//下界
int *t_U;//辅助数组，用来记录BB里边每次递归里的U

int **block;//用来记录已经开辟的内存块，每一块大小设定为2G
int cur_block = -1;//当前内存块编号
int block_size = 0;//已开辟内存块数量

unordered_set<unsigned long long> edge_hash;//边哈希，用来快速判断顶点之间是否有边，避免使用邻接表结构

void freeAll();//释放内存函数

bool cmp (Vertex_sort v1, Vertex_sort v2)//用于顶点排序，从大到小排序
{
    return v1.degree > v2.degree;
}

double get_utime()//获取已经运行时间，从程序运行开始计时
{
    struct rusage utime;
    getrusage(RUSAGE_SELF, &utime);
    return (double) (utime.ru_utime.tv_sec + (double)utime.ru_utime.tv_usec / 1000000);
}

unsigned long long encode_pairID(unsigned long long  v1, unsigned long long  v2)//用来编码边的哈希值
{
    unsigned long long n1,n2;
    unsigned long long pairID;
    if(v1 < v2)
    {
        n1 = v1; n2 = v2;
    }
    else
    {
        n1 = v2; n2 = v1;
    }
    pairID = ((n1 + n2 + 1) * (n1 + n2) >> 1) + n2;
    return pairID;
}

int edge_is(int m, int n)//判断m,n间是否有边，有返1，否0
{
    unsigned long long id = encode_pairID(m, n);
    if(edge_hash.count(id))
        return 1;
    else
        return 0;
}

int* expand_memory()//扩展内存块，如果无法拓展内存，退出程序，否则返回下一块内存的首地址
{
    if(cur_block == block_size - 1)
    {
        block[block_size] = (int *)malloc(mem_2G * sizeof(int));
        if(block[block_size] == NULL)
        {
            freeAll();
            exit(0);
        }
        block_size++;
    }
    cur_block++;
    return  block[cur_block];
}

void readGraph(char* File_Name)//读入原始图，存入相应的数据结构
{

    ifstream FIC;
    FIC.open(File_Name);
    if(FIC.fail())
    {
        printf("### Error open, File_Name %s\n", File_Name);
    }
    char strReading[1024];
    char p[10], temp[10];
    FIC >> p >> temp >> v_n >> e_n;

    neighbor = (int **)malloc(v_n * sizeof(int**));
    block = (int **)malloc(100 * sizeof(int**));
    neighbor_len = (int *)malloc(v_n * sizeof(int));
    vertex_sort= (Vertex_sort *)malloc(v_n *sizeof(Vertex_sort));
    edge = (Edge *)malloc(e_n * sizeof(Edge));
    temp_array = (int *)malloc(v_n * sizeof(int));
    temp_mark = (int *)malloc(v_n * sizeof(int));
    temp_index = (int *)malloc(v_n * sizeof(int));
    vertex = (Vertex *)malloc(v_n * sizeof(Vertex));
    neighbor_in_solution = (int *)malloc(v_n * sizeof(int));
    t_U = expand_memory();
    best_solution = (int *)malloc(v_n *sizeof(int));

    remaining_vertex = new Array(v_n);///////////////////////////////////////////
    crusol = new Array(v_n);
    U = new Array(v_n);

    for(int i = 0; i < v_n; i++)
    {
        vertex[i].degree = 0;
        vertex[i].state = State::Candidate;
        vertex[i].num_in_c = 0;
        vertex[i].is_in_search = 0;
        vertex[i].score = 0;
    }

    int a,b;
    int count = 0;
    unsigned long long id;

    while(FIC >> temp >> a >> b)
    {
        if(strcmp(temp, "e") ==  0)
        {
            a--;
            b--;
            vertex[a].degree++;
            vertex[b].degree++;
            edge[count].v1 = a;
            edge[count].v2 = b;
            id = encode_pairID(a, b);
            edge_hash.insert(id);
            count++;
        }
    }

    for(int i = 0; i < v_n; i++)
    {
        if(vertex[i].degree > max_degree)
            max_degree = vertex[i].degree;
        vertex[i].neighbor = (int *)malloc(vertex[i].degree * sizeof(int));
        neighbor[i] = (int *)malloc(vertex[i].degree * sizeof(int));
        remaining_vertex->insert_element(i);
    }

    for(int i = 0; i < v_n; i++)
        vertex[i].degree = 0;
    for(int e  = 0; e < e_n; e++)
    {
        a = edge[e].v1;
        b = edge[e].v2;
        vertex[a].neighbor[vertex[a].degree] = b;
        vertex[b].neighbor[vertex[b].degree] = a;
        vertex[a].degree++;
        vertex[b].degree++;
    }
    //计算顶点score值
    for (size_t i = 0; i < v_n; i++)
    {
        vertex[i].score = vertex[i].degree + 1;//自身也要添加
    }
    remain_num = v_n;
    uncover_num = v_n;
    free(edge);
    FIC.close();
    printf("### the number of vertex:\t%d\n### the number of edge:\t\t%d\n### maximum degree:\t\t%d\n### time limit: \t\t%ds\n",
           v_n, e_n, max_degree, run_time);
}

void freeAll()//释放内存函数
{
    free(neighbor_len);
    free(temp_array);
    free(temp_mark);
    free(temp_index);
    free(neighbor_in_solution);
    free(best_solution);
    for(int i = 0; i <100; i++)
        free(block[i]);
    free(block);
    for(int i = 0; i < v_n; i++)
    {
        free(vertex[i].neighbor);
        free(neighbor[i]);
    }
    free(vertex);
    free(neighbor);
    free(vertex_sort);
}

void update_best_solution()//更新最优解
{
    best_solution_size = crusol->size();
    for(int i = crusol->begin(); i < crusol->size(); i++)
        best_solution[i] = crusol->element_at(i);
    LB = best_solution_size;
}

bool check_solution()//检验最优解是否可行，是返回1，否则返0
{
    int i, j, adj_num, v, u;
    for(i = 0; i < best_solution_size; i++)
    {
        adj_num = 0;
        v = best_solution[i];
        for(j = 0; j < best_solution_size; j++)
        {
            u = best_solution[j];
            if(v ==u)
                continue;
            if(edge_is(v,u))
                adj_num++;
        }
        if(adj_num < best_solution_size - para_k)
            return false;
    }
    return true;
}

void printf_solution()//打印最优解，这边只输出当前最优解的规模
{
    if(check_solution())
    {
        printf("after checking , the solution is correct,  solution size: %d, time: %f\n", best_solution_size, get_utime());
    }
    else
    {
        printf("the solution found is wrong\n");
    }

}

void reduce_v(int v)//如果v的score为0,reduce
{
    if (vertex[v].score == 0)
        reduce(v);
}

void fix_vertex(int v)//固定顶点v到最优解
{
    if (vertex[v].state == State::Fixed)
        return;
    fix_num++;
    fix(v);
    vertex[v].score = -vertex[v].score;
    if (vertex[v].num_in_c == 0)
        uncover_num--;
    for (size_t i = 0; i < neighbor_len[v]; i++)
    {
        int v_neighbor = vertex[v].neighbor[i];
        if (vertex[v_neighbor].num_in_c == 0)
            uncover_num--;
        vertex[v_neighbor].num_in_c++;
        if (vertex[v].num_in_c == 0)
        {
            //v首次被支配，修改其邻居顶点
            vertex[v_neighbor].score -= 1;
            reduce_v(v_neighbor);
        }
        else if (vertex[v].num_in_c == 1 && vertex[v_neighbor].state == State::Fixed)
        {
            //v_neighbor存在在候选解，score值为负
            vertex[v_neighbor].score += 1;
            reduce_v(v_neighbor);
        }
        if (vertex[v_neighbor].state == State::Fixed)
        {
            if (vertex[v_neighbor].num_in_c == 2)
            {
                //v_neighbor已经在候选解，v加入后，v_neighbor又被支配了一次
                vertex[v_neighbor].score += 1;
                reduce_v(v_neighbor);
            }
            continue;
        }
        //处理二层邻居
        int cnt = 0;
        int s = 0;
        for (size_t j = 0; j < neighbor_len[v_neighbor]; j++)
        {
            int v_n_n = vertex[v_neighbor].neighbor[j];//二层邻居
            if (v_n_n == v)
            {
                continue;
            }
            if (vertex[v_n_n].state == State::Fixed)
            {
                s = v_n_n;
                cnt++;
            }
        }
        if (vertex[v_neighbor].state == State::Fixed)
        {
            s = v_neighbor;
            cnt++;
        }
        if (cnt == 0)
        {
            //v是闭邻居里第一个加入候选解的
            vertex[v_neighbor].score -= 1;
            reduce_v(v_neighbor);
            
        }
    }
    
}

void reduce_graph_1()//reduce边缘顶点
{
    for (size_t i = 0; i < v_n; i++)
    {
        if (vertex[i].degree == 1)
        {
            int v_neighbor = vertex[i].neighbor[0];
            if (vertex[i].state != State::Fixed)
            {
                reduce(i);
            }
        }
    }
}

/*void reduce_edge_from_neighbor(int v, int mark)//update the neighbor list//更新边邻居，因为复杂度有些偏高，所以优化的时候，把这部分改掉了
{
    int i,j;
    int u;
    if(mark == 1)
    {
        for(i = 0; i < vertex[v].degree; i++)
        {
            u = vertex[v].neighbor[i];
            for(j = 0; j < vertex[u].degree; j++)
                if(vertex[u].neighbor[j] == v)
                    break;
            vertex[u].degree--;
            vertex[u].neighbor[j] = vertex[u].neighbor[vertex[u].degree];
        }
        vertex[v].degree = 0;
    }
    else if(mark == 2)
    {
        for(i = 0; i < neighbor_len[v]; i++)
        {
            u = neighbor[v][i];
            for(j = 0; j < neighbor_len[u]; j++)
                if(neighbor[u][j] == v)
                    break;
            neighbor_len[u]--;
            neighbor[u][j] = neighbor[u][neighbor_len[u]];
        }
        //neighbor_len[v] = -1;
    }
}*/

/*void reduce_graph_1()//初始化里边的reduce1，用的是团那边类似的初始化方法求下界
{
    int *degree_counter, *where, *candidate;  
    int node, p1, i, j, h, k, t, neighbors,_tsize;
    int u,v;
    int cur_degree = 0;
    degree_counter = temp_mark;
    where = temp_index;
    candidate = temp_array;

    for(i = 0; i <= max_degree; i++)
        degree_counter[i] = 0;
    for(node = remaining_vertex->begin(); node < remaining_vertex->size(); node++)
    {
        v = remaining_vertex->element_at(node);
        degree_counter[vertex[v].degree]++;
    }
    j = 0;
    for(i = 0; i <= max_degree; i++)//degree_counter每一个度大小的起始位置指针
    {
        k = degree_counter[i];
        degree_counter[i] = j;
        j += k;
    }

    for(node = remaining_vertex->begin(); node < remaining_vertex->size(); node++)
    {
        v = remaining_vertex->element_at(node);
        t = degree_counter[vertex[v].degree];
        degree_counter[vertex[v].degree]++;
        candidate[t] = v;
        where[v]  = t;
    }

    for(i = max_degree; i > 0; i--)
        degree_counter[i] = degree_counter[i - 1];

    degree_counter[0] = 0;
    p1 = 0;
    cur_degree = vertex[candidate[0]].degree;

    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)
    {
        v = remaining_vertex->element_at(i);
        neighbor_len[v] = vertex[v].degree;
    }

    while(p1 < remaining_vertex->size())
    {
        node = candidate[p1];
        if(p1 < remaining_vertex->size() - 1  && neighbor_len[node] == neighbor_len[candidate[p1+1]]) //如果数组后面的一个度数和他相等，更新一下该度数大小的指针
            degree_counter[neighbor_len[node]] = p1 + 1;
        if(neighbor_len[node] + para_k - 1 >= remaining_vertex->size() - p1 - 1)   //团：该顶点的度数==后面的顶点度数之和  则该顶点和后面顶点构成一个团
        {
            crusol->clear();
            for(i = p1; i < remaining_vertex->size(); i++)
                crusol->insert_element(candidate[i]);
            break;
        }

        for(i = 0; i < vertex[node].degree; i++)   //删除node，之后更新他的所有邻居
        {
            neighbors = vertex[node].neighbor[i];  //	原图邻居   
            if(where[neighbors] > p1)   //只需要判断原图邻居 在P1后面的位置上的顶点
            {
                t = where[neighbors];   //neighbors的地址为t
                h = degree_counter[neighbor_len[neighbors]];   //neighbors的度数的首地址为h
                k = candidate[h]; //neighbors的度数的首地址中的变量为k

                candidate[h] = neighbors;  //neighbors和k互换位置
                where[neighbors] = h;
                candidate[t] = k;
                where[k]= t;

                degree_counter[neighbor_len[neighbors]]++;
                neighbor_len[neighbors]--;
                if(neighbor_len[neighbors] != neighbor_len[candidate[h-1]])// 若度数地址 存在 3 3 3 5  这种 中间没有度数4的情况，需要更新
                    degree_counter[neighbor_len[neighbors]] = h;
            }
        }
        p1++;
    }

    if(crusol->size() > best_solution_size)
        update_best_solution();

    queue<int> remove_que;
    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)
    {
        v = remaining_vertex->element_at(i);
        neighbor_len[v] = vertex[v].degree;
        if(vertex[v].degree + para_k <= best_solution_size ||  remaining_vertex->size() <= best_solution_size)
        {
            vertex[v].state = State::Delete;
            remove_que.push(v);
        }
    }

    while(!remove_que.empty())//迭代删除顶点
    {
            v = remove_que.front();
            remove_que.pop();
            remaining_vertex->delete_element(v);
            for(i = 0; i < vertex[v].degree; i++)
            {
                u = vertex[v].neighbor[i];
                neighbor_len[u]--;
                if(!vertex[u].state && (neighbor_len[u] + para_k <= best_solution_size ||  remaining_vertex->size() <= best_solution_size))
                {
                    vertex[u].state = -1;
                    remove_que.push(u);
                }
            }
    }

    int mark = 1;
    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)//更新剩下点的邻居信息
    {
        v = remaining_vertex->element_at(i);
        mark = 1;
        for(j = 0; j < vertex[v].degree; j++)
        {
            u = vertex[v].neighbor[j];
            if(vertex[u].state)
            {
                mark = 0;
                vertex[v].degree--;
                vertex[v].neighbor[j] = vertex[v].neighbor[vertex[v].degree];
                j--;
            }
        }
        if(!mark && U->is_in_array(v))//U用来存放在初始化过程中顶点度没有变化的顶点集合，如果在reduce过程中U里顶点度发生改变，则从U中删除
            U->delete_element(v);
    }
}*/


/*void  reduce_graph_2(Array *Uset)//用在初始化里边的reduce2，Uset=remaining_vertex
{
    int i,v,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int ssize;
    ssize = 0;//初始化时当前解大小为0
    queue<int> remove_que;
    for(i = Uset->begin(); i < Uset->size(); i++)
    {
        u1 = Uset->element_at(i);
        neighbor_len[u1] = 0;
    }

    for(i = Uset->begin(); i < Uset->size(); i++)
    {
        v = Uset->element_at(i);//v是这次迭代中要判断是否移除的顶点
        if(U->is_in_array(v))//U用来存放在初始化过程中顶点度没有变化的顶点集合，如果v在上次reduce1，reduce2过程中度没有发生改变，则跳过这个顶点
            continue;
        kv = para_k - 1;//初始化时，kv值都是k-1
        rv = LB - kv - ssize;
        rd = vertex[v].degree;
        neighbor_len[v] = 0;
        for(j = 0; j < rd; j++)//初始化G'里边的邻居信息
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))
                continue;
            neighbor_len[u1] = 1;
            neighbor[u1][0] = v;
            //neighbor_len[v]++;
            neighbor[v][neighbor_len[v]++] = u1;
        }
        if(neighbor_len[v]  < rv)//如果顶点v的度小于rv，则删除这个顶点
        {
            Uset->delete_element(v);
            //reduce_edge_from_neighbor(v, 1);
            i--;
            continue;
        }
        while(!remove_que.empty())//清空移除队列
            remove_que.pop();
        for(j = 0; j < rd; j++)//构建子图G'
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))//如果点v已经被移除了，则跳过
                continue;
            for(k = j + 1; k < rd; k++)
            {
                u2 = vertex[v].neighbor[k];
                if(!Uset->is_in_array(u2))
                    continue;
                if(edge_is(u1, u2))
                {
                    neighbor[u1][neighbor_len[u1]++] = u2;
                    neighbor[u2][neighbor_len[u2]++] = u1;
                }
            }
            ku = para_k - 1;
            cu = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
            temp_index[u1] = neighbor_len[u1];//这边用temp_index来记录G'中每个点在G'中的度，方便后边迭代时更新
            if(cu < rv)//如果cu<rv，则加入移除队列
            {
                temp_index[u1] = 0;//temp_index=0表示这个点已经移除出G'
                remove_que.push(u1);
            }
        }

        temp_index[v] = neighbor_len[v];//顶点v的temp_index值等于它在G'中的度
        while(!remove_que.empty())//迭代删除G'中的顶点
        {
            u1 = remove_que.front();
            remove_que.pop();
            for(j = 0; j < neighbor_len[u1]; j++)
            {
                u2 = neighbor[u1][j];
                if(!temp_index[u2])//如果顶点u2已经移除G'，则跳过
                    continue;
                temp_index[u2]--;//更新邻居顶点的度数
                if(temp_index[u2] + ku < rv && u2 != v)//如果cu<rv并且u2不是当前的v
                {
                    temp_index[u2] = 0;//temp_index=0表示这个点已经移除出G'
                    remove_que.push(u2);
                }
            }
            if(temp_index[v] < rv)//如果顶点v的度已经小于rv，则删除v，同时跳出迭代
            {
                Uset->delete_element(v);
                i--;
                break;
            }
        }
    }

    U->clear();//清空U
    int mark = 1;
    for(i = Uset->begin(); i < Uset->size(); i++)//更新未删除顶点的邻居信息
    {
        v = Uset->element_at(i);
        mark = 1;
        for(j = 0; j < vertex[v].degree; j++)
        {
            u1 = vertex[v].neighbor[j];
            if(!Uset->is_in_array(u1))
            {
                mark = 0;
                vertex[v].degree--;
                vertex[v].neighbor[j] = vertex[v].neighbor[vertex[v].degree];
                j--;
            }
        }
        if(mark)//如果U在reduce2过程中度数没有改变，则加入U
            U->insert_element(v);
    }
}*/

void  reduce_graph_in_BB(int* &begin, int* &end, int rev)//BB过程中的reduce2，begin是此次递归过程中U的起始地址，end结束地址，remnum是此次reduce过程中reduce1加reduceS约减掉的顶点
{
    int i,v,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int mark = 1;
    int ssize;
    int *ii;
    int t_LB = -1;
    queue<int> remove_que;
    ssize = S->size();//当前解大小

    for(ii = begin; ii < end; ii++)
        temp_index[*ii] = 1;//用来判断是否在begin---end范围内，即判断顶点v是否在U中，其实这语句不加也行，因为temp_index值肯定是1

    for(i = 0; i < vertex[rev].degree; i++)
    {
        u1 = vertex[rev].neighbor[i];
        if(temp_index[u1])//在begin，end范围内，且是rev的邻居。reduce2只针对rev的邻居进行判断
            temp_index[u1] = 2;
    }

    for(ii = begin; ii < end; ii++)
    {
        mark = 1;
        v  = *ii;
        if(temp_index[v] < 2)
            continue;
        kv = para_k - 1 - (ssize - neighbor_in_solution[v]);//顶点v的kv值   kv是可以与v没有边的顶点个数
        /*if(kv < 0)//如果kv<0，则直接移除这个顶点
        {
            cout << "@";
            delete_i(ii, end);
            ii--;
            temp_index[v] = 0;//移除点以后，temp_index值更新为0
            /*for(i = 0; i < vertex[v].degree; i++)
            {
                u1 = vertex[v].neighbor[i];
                if(temp_index[u1])
                    vertex[u1].state--;
            }
            continue;
        }*/
        rv = LB - kv - ssize;//顶点v的rv值  rv是在候选顶点中应该与v有边顶点的个数
        rd = vertex[v].degree;
        neighbor_len[v] = 0;
        for(j = 0; j < rd; j++)//初始化G‘的信息
        {
            u1 = vertex[v].neighbor[j];
            if(!temp_index[u1])
                continue;
            neighbor[u1][0] = v;
            neighbor_len[u1] = 1;
            neighbor[v][neighbor_len[v]++] = u1;
        }

        if(neighbor_len[v]  < rv)//如果顶点v的度小于rv，直接移除
        {
            delete_i(ii, end);
            ii--;
            temp_index[v] = 0;
          /*  for(i  = 0; i < neighbor_len[v]; i++)
            {
                u1 = neighbor[v][i];
                vertex[u1].state--;
            }*/
            continue;
        }

        while(!remove_que.empty())
            remove_que.pop();

        for(j = 0; j < rd; j++)//构建子图G'
        {
            u1 = vertex[v].neighbor[j];
            if(!temp_index[u1])
                continue;
            for(k = j + 1; k < rd; k++)
            {
                u2 = vertex[v].neighbor[k];
                if(!temp_index[u2])
                    continue;
                if(edge_is(u1, u2))
                {
                    neighbor[u1][neighbor_len[u1]++] = u2;
                    neighbor[u2][neighbor_len[u2]++] = u1;
                }
            }

            ku = para_k - 1 - (ssize - neighbor_in_solution[u1]);  //ku是u1还可以与多少个顶点没有边
            temp_array[u1] = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
            //temp_array[u1] = neighbor_len[u1] + ku;
            if(temp_array[u1] < rv )
            {
                temp_array[u1] = 0;//temp_arry=0说明这个点已经移除出G'
                remove_que.push(u1);
            }
        }

        temp_array[v] = neighbor_len[v];//顶点v的temp_array值为这个点的度数
        while(!remove_que.empty())//迭代更新G'
        {
            u1 = remove_que.front();
            remove_que.pop();
            temp_array[v]--;//v的邻居数--
            for(j = 1; j < neighbor_len[u1]; j++)//邻居中从1下标1开始，因为下标0存的是v
            {
                u2 = neighbor[u1][j];
                if(!temp_array[u2])//如果u2已经移除出G'，则跳过
                    continue;
                temp_array[u2]--;//更新顶点u2的cu值
                if(temp_array[u2] > temp_array[v])//如果cu值大于v的度，因为cu=min（du+ku，G'中的顶点数）
                    temp_array[u2] = temp_array[v];
                if(temp_array[u2] < rv)//如果顶点u2的cu值小于rv，移除u2
                {
                    temp_array[u2] =0;
                    remove_que.push(u2);
                }
            }
            if(temp_array[v] < rv)//如果顶点v的度小于rv，从U中删除v
            {
                delete_i(ii, end);
                ii--;//因为删除操作是把最后一个点覆盖到当前位置，所以下标指针要前移1位
                temp_index[v] = 0;//用来表示v是否在U中，不在=0
            /*    for(i = 0; i < neighbor_len[v]; i++)
                {
                    u1 = neighbor[v][i];
                    vertex[u1].state--;
                }*/
                //mark = 0;
                break;
            }
        }
    }

    for(ii = begin; ii < end; ii++)//计算U中顶点在U中的度，存入state变量中，因为reduce2过程中没有动态更新度，动态更新更费时
    {
        v = *ii;
        vertex[v].state = 0;
        for(i = 0; i < vertex[v].degree; i++)
        {
            if(temp_index[vertex[v].neighbor[i]])
                vertex[v].state++;
        }
    }
}

void reduce_graph()//初始化约减函数
{
    int rn = remaining_vertex->size();
    printf("best solution\tcurrent solution\tremaining vertex\n");
    printf("%8d\t%8d\t\t%8d\n", 0, 0, v_n);
    while(true)//迭代reduce1和reduce2，直至无法约减
    {
        reduce_graph_1();
        printf("%8d\t%8d\t\t%8d(reduce1)\n", best_solution_size, crusol->size(), remaining_vertex->size());
        reduce_graph_2(remaining_vertex);
        printf("%8d\t%8d\t\t%8d(reduce2)\n", best_solution_size, crusol->size(), remaining_vertex->size());

        if(remaining_vertex->size() < rn)
            rn = remaining_vertex->size();
        else
            break;
    }

}

int calculate_upbound(Array *ary)//计算初始化约减以后的上界，其实没有用到，可以注释掉。ary是要计算上界的数据集，数据集格式是Array
{

    int i,vd,v,j;
    for(i = ary->begin(); i < ary->size(); i++)
    {
        vd = 0;
        v = ary->element_at(i);
        vertex_sort[i].index = v;
        vertex_sort[i].degree = vertex[v].degree;
    }
    sort(vertex_sort, vertex_sort + ary->size(), cmp);
    for(i = ary->begin(); i < ary->size(); i++)
    {
        vd = vertex_sort[i].degree;
        if(i + 1 >= vd + para_k - 1)
            break;
    }
    return i + 1;
}

int calculate_upbound_in_BB(int *begin, int *end)//BB过程中计算上界，begin是此次递归中U的起始地址，end是终止地址
{
    int i = 0,vd,v,j;
    int *ii;
    int tt = 0;
    for(ii = begin; ii < end; ii++)
    {
        v = *ii;
        //vertex_sort[i].index = v;
        //vertex_sort[i].degree = vertex[v].state + para_k  - (S->size() - neighbor_in_solution[v]);//顶点v在U中的上界是du+kv+1，+1是因为需要算自身
        //vertex_sort[i].degree = vd + para_k  - (S->size() - neighbor_in_solution[v]);
        neighbor_len[v] = vertex[v].state;
        if(neighbor_len[v] + para_k + neighbor_in_solution[v] > LB)
            tt++;
        temp_index[v] = 0;
    }
    return tt;

    //sort(vertex_sort, vertex_sort + i, cmp);//排序，寻找上界，即找d，使得集合里边kplex度大于等于d的顶点数量大于等于d
    //for(ii = begin; ii < end; ii++)
        //temp_index[*ii] = 0;
   // for(i = 0; i < remaining_vertex->size(); i++)
       // if(temp_index[remaining_vertex->element_at(i)])
          //  cout << "$";
    /*for(j = 0; j < i; j++)
    {
        vd = vertex_sort[j].degree;
        if(j + 1 >= vd)
            break;
    }*/

   // if(j + 1 < LB && end-begin > LB)
    //cout << j+1 << "&" << end - begin << "&" << tt <<"&" << LB - S->size() << " ";
    //return j + 1;
    //return tt;

}

void reduceU1(int* &tbegin, int* &tend)
{
    int *ii, v, i, u;
    queue<int> remove_que;
    for(ii = tbegin; ii < tend; ii++)//标记begig，end中的顶点，如果在begin-end过程中，则temp_index值为1
        temp_index[*ii] = 1;

    for(ii = tbegin; ii < tend; ii++)
    {
        v = *ii;
        vertex[v].state = 0;
        for(i = 0; i < vertex[v].degree; i++)
            if(temp_index[vertex[v].neighbor[i]] > 0)
                vertex[v].state++;//计算U中的顶点在U中度
        if(vertex[v].state + neighbor_in_solution[v] + para_k <= LB)//如果顶点v针对LB的上界展望小于LB，则删除这个顶点
        {//候选顶点中v的邻居+S解中v的邻居+k<=LB
            remove_que.push(v);
            temp_index[v] = 2;
        }
    }

    while(!remove_que.empty())//迭代删除U中的顶点
    {
        v = remove_que.front();
        remove_que.pop();
        temp_index[v] = 0;
        for(i = 0; i < vertex[v].degree; i++)
        {
            u = vertex[v].neighbor[i];
            if(temp_index[u] != 1)
                continue;
            vertex[u].state--;
            if(vertex[u].state + neighbor_in_solution[u] + para_k <= LB)
            {
                temp_index[u] = 0;
                remove_que.push(u);
            }
        }
    }

    for(ii = tbegin; ii < tend; ii++)//通过temp_index值更新U，因为U删除的时候是根据地址删除的，所以需要通过遍历来做
    {
        v = *ii;
        if(!temp_index[v])
        {
            delete_i(ii, tend);
            ii--;
        }
    }
}

void reduceV(int* &begin, int* &end, int v)//begin，end用来表示U的首尾指针，v是此次加入或删除解集的顶点
{
    int i,j,k,p,u1,u2,tnode;
    int rv, rd, cu, kv, ku;
    int ssize;
    int *ii;
    queue<int> remove_que;
    ssize = S->size();
    kv = para_k - (ssize - neighbor_in_solution[v]);//顶点v的kv值
    rv = LB - kv - ssize;//顶点v的rv值
    rd = vertex[v].degree;
    neighbor_len[v] = 0;
    for(j = 0; j < rd; j++)//初始化邻居信息
    {
        u1 = vertex[v].neighbor[j];
        if(!temp_index[u1])
            continue;
        neighbor[u1][0] = v;
        neighbor_len[u1] = 1;
        neighbor[v][neighbor_len[v]++] = u1;
    }
    for(j = 0; j < rd; j++)//构建v的诱导子图G'
    {
        u1 = vertex[v].neighbor[j];
        if(!temp_index[u1])
            continue;
        for(k = j + 1; k < rd; k++)
        {
            u2 = vertex[v].neighbor[k];
            if(!temp_index[u2])
                continue;
            if(edge_is(u1, u2))
            {
                neighbor[u1][neighbor_len[u1]++] = u2;
                neighbor[u2][neighbor_len[u2]++] = u1;
            }
        }
        ku = para_k - 1 - (ssize - neighbor_in_solution[u1]);
        temp_array[u1] = neighbor_len[u1] + ku < neighbor_len[v]? neighbor_len[u1] + ku : neighbor_len[v];
        if(temp_array[u1] < rv)//如果G'中u的cu值小于rv，则删除点u
        {
            temp_array[u1] = 0;
            remove_que.push(u1);
            //cout << "%";
        }
    }
    temp_array[v] = neighbor_len[v];

    while(!remove_que.empty())//迭代删除
    {
        u1 = remove_que.front();
        remove_que.pop();
        temp_array[v]--;//用来表示是否在G'中，如果在的话表示顶点cu值，只有顶点v表示的它的度
        temp_index[u1] = 0;
        for(j = 0; j < vertex[u1].degree; j++)
        {
            u2 = vertex[u1].neighbor[j];//更新u1邻居的在U中的度，不能用neighbor_len信息，因为neighbor_len是在G’的信息，不是U中信息
            if(temp_index[u2])
                vertex[u2].state--;
        }
        for(j = 1; j < neighbor_len[u1]; j++)//更新u1邻居的cu值
        {
            u2 = neighbor[u1][j];
            //vertex[u2].state--;
            if(!temp_array[u2])
                continue;
            temp_array[u2]--;
            if(temp_array[u2] > temp_array[v])
                temp_array[u2] = temp_array[v];
            if(temp_array[u2] < rv)
            {
                temp_array[u2] = 0;
                remove_que.push(u2);
            }
        }
    }

    for(ii = begin; ii < end; ii++)//通过temp_index来更新U
    {
        u1 = *ii;
        if(!temp_index[u1])
        {
            delete_i(ii, end);
            ii--;
        }
    }
}

int is_S;
int nn= 0;
void update_candidate(int* &tbegin, int* &tend, int mark, int rev)//更新约减U，tbegin是此次U的起始位置，tend是终止位置，其实相当于上边函数的begin和end值，mark是push/pop的指示器1表示是pop，0表示push, rev表示此次加入或删除的顶点
{
    int tt1,tt2 = 0,tt3 = 0,tt4;
    int v, i, j;
    int *ii;
    //tt1 = tend-tbegin;
    if(mark)//pop过程，进行reduce1和reduce2
    {
        //cout << "\npop:";
        reduceU1(tbegin, tend);
        //tt2 = tend-tbegin;
        reduce_graph_in_BB(tbegin, tend, rev);
    }
    else//push过程，进行reduceS，reduce1，reduceV过程
    {
        //cout << "\npush:";
        for(ii = tbegin; ii < tend; ii++)
            temp_mark[*ii] = 0;
        is_S = 0;
        for(i = S->begin(); i < S->size(); i++)//遍历S，判断是否有点已经到达kplex的临界条件，即kv==0
        {
            v = S->element_at(i);
            if(neighbor_in_solution[v] + para_k == S->size())  //不能再选v的非邻居顶点了
            {
                for(j = 0; j < vertex[v].degree; j++)
                        temp_mark[vertex[v].neighbor[j]]++;
                is_S++;
            }
        }
        for(ii = tbegin; ii < tend; ii++)
        {
            v = *ii;
            if((is_S && temp_mark[v] != is_S) || neighbor_in_solution[v] + para_k <= S->size())
            {//temp_mark[v] != is_S  这个v与临界顶点都有边，可以不删
                delete_i(ii, tend);
                //temp_index[v] = 0;
                ii--;
            }
        }

        //tt2 = tend - tbegin;
        reduceU1(tbegin, tend);
       //tt3 = tend - tbegin;
        reduceV(tbegin, tend, rev);
    }
    /*if(mark)
    {
        reduce_graph_in_BB(tbegin, tend, rev);
    }
    else
    {
        reduceV(tbegin, tend, rev);
    }*/
    //tt4 = tend - tbegin;
    //if(mark == 0 && tt2 == tt1 && tt3 < tt2)
        //cout << tt1 << "@" << tt2 << "@" << tt3 << "@" << tt4 << " ";
}

int ct =  0;//ct用来记录递归函数调用次数，即递归展开的节点数

int BB(int* &tbegin, int* &tend, int *tv, int nt)//分支界限过程，tbegin，tend分别是此次递归中U的起始和终止地址，vt是U中度最大的顶点地址，nt是递归深度
{

    if(get_utime() > run_time)//如果超出time_limit，则释放内存后退出程序
    {
        freeAll();
        printf("out of time\n");
        exit(0);
    }
    int u, i, upper, lower, maxd, v, vd, d;
    int *maxv;
    int *ii;
    int *begin, *end;
    int *cur_ind = block[cur_block];//记录当前内存块编号
    int t_cur_block = cur_block;
    if(tbegin < tend)
    {
        ct++;
        nt++;
        v = *tv;
        delete_i(tv, tend);//从U中删除v      begin-end里面是所有候选顶点     S中存的当前解中顶点
        S->insert_element(v);//v加入S
        for(i = 0; i < vertex[v].degree; i++)//更新顶点与S的相邻数组
            neighbor_in_solution[vertex[v].neighbor[i]]++;  //当前解中与该顶点是邻居的个数
        begin = tend + 1;//开辟下一次递归的U，可以认为是t_U    //begin保存的是递归下次存放的候选顶点的开始位置
        if(tend - tbegin > cur_ind + mem_2G - tend - 2)//如果当前内存块已经无法开辟U，则开辟一块新的内存块
        {
			//tend - tbegin当前用的多少； cur_ind + mem_2G - tend剩余还有多少
            begin = expand_memory();
        }
        end = begin;
        for(ii = tbegin; ii < tend; ii++)//将这次的U拷贝到下一次的U（t_U）
        {
            insert_v(end, *ii);  // 将候选顶点复制到begin到end中
        }

       // cout << end - begin << "&";
        //reduceV(begin, end, v);
        //cout << end - begin << " ";
        update_candidate(begin, end, 0, v);//更新约减t_U，传入的参数0表示是push过程，v表示此次加入解集的顶点
        upper = calculate_upbound_in_BB(begin, end) + S->size();//计算t_U的上界，需要最后加上S的大小
        if(upper > LB)
        {
            maxd = -1;
 	    int minkv=para_k;
            for(ii = begin; ii < end; ii++)//寻找t_U中度最大的顶点
            {
                u = *ii;
                //d = neighbor_len[u] + para_k - 1 - (S->size() - neighbor_in_solution[u]);//顶点在t_u的度+kv
                d = neighbor_len[u];//顶点在t_U中的度
                int kkk=para_k - 1 - (S->size() - neighbor_in_solution[u]);
		if(kkk==0 ){
		   if(minkv==0 && d>maxd || minkv>0){
		   maxv = ii;
		   minkv= kkk;
		   maxd = d;
		     }

		}

		else if(minkv>0 && d > maxd)
                {

                    maxd = d;
                    maxv = ii;
		   minkv=kkk;
                }

/*if(d > maxd){
  maxd = d;
                    maxv = ii;

}*/

            }

            lower = BB(begin, end, maxv, nt);
            if(lower > LB)//如果找到更好的解，更新解
            {
                LB = lower;
                update_best_solution();
                printf_solution();
            }
        }
        S->delete_element(v);//从S中删除v

        for(i = 0; i < vertex[v].degree; i++)//更新顶点与S的相邻数组
            neighbor_in_solution[vertex[v].neighbor[i]]--;
        cur_block = t_cur_block;//恢复当前内存块
        begin = tend + 1;
        if(tend - tbegin > cur_ind + mem_2G - tend - 2)
        {
            begin = expand_memory();
        }
        end = begin;
        for(ii = tbegin; ii < tend; ii++)//将这次的U拷贝到下一次的U（t_U）
            insert_v(end, *ii);
        update_candidate(begin, end, 1, v);//更新约减t_U，传入的参数1表示是pop过程，v表示此次从解集里删除的顶点
        upper = calculate_upbound_in_BB(begin, end) + S->size();
        if(upper > LB)
        {
            maxd = -1;
	    int minkv=para_k;
            for(ii = begin; ii < end; ii++)//寻找t_U中度最大的顶点
            {
                u = *ii;
                //d = neighbor_len[u] + para_k - 1 - (S->size() - neighbor_in_solution[u]);//顶点在t_u的度+kv
                d = neighbor_len[u];//顶点在t_U中的度数
                int kkk=para_k - 1 - (S->size() - neighbor_in_solution[u]);
		if(kkk==0 ){
		   if(minkv==0 && d>maxd || minkv>0){
		   maxv = ii;
		   minkv= kkk;
		   maxd = d;
		     }

		}

		else if(minkv>0 && d > maxd)
                {

                    maxd = d;
                    maxv = ii;
		   minkv=kkk;
                }

/*if(d > maxd){
  maxd = d;
                    maxv = ii;

}*/

            }
            lower = BB(begin, end, maxv, nt);
            if(lower > LB)
            {
                LB = lower;
                update_best_solution();
                printf_solution();
            }
        }
        return LB;
    }
    else
    {
        //cout << nt << " ";
        return S->size();
    }

}

void search()
{
    int i, v,  maxd = -1;
    int *begin, *end;
    int *maxv;
    end = t_U;
    begin = t_U;
    int en = 0;
    for(i = remaining_vertex->begin(); i < remaining_vertex->size(); i++)
    {
        v = remaining_vertex->element_at(i);
        //U->insert_element(v);
        neighbor_in_solution[v] = 0;
        insert_v(end, v);//初始化BB第一次的U
        temp_index[v] = 0;
        en += vertex[v].degree;
        if(vertex[v].degree > maxd)//记录度最大的顶点
        {
            maxd = vertex[v].degree;
            maxv = end;
        }
    }
    crusol->clear();//清空当前解S
    S = crusol;
   // cout << "after reduce, the edge num = " << en / 2 << "\n";
    LB = BB(begin, end, maxv, 0);
}

int main(int argc, char *argv[])
{
    char* File_Name;
    File_Name = argv[1];//读入文件名
    para_k = atoi(argv[2]);//参数k
    run_time = atoi(argv[3]);//time_limit，运行时间阈值
    printf("%s: \n", File_Name);
    readGraph(File_Name);//读入原图
    read_time = get_utime();
    printf("----------------------------readGraph finish----------------------------\n");
    reduce_graph();//约减原图
    init_time = get_utime() - read_time;
    LB = best_solution_size;
    //UB = calculate_upbound(remaining_vertex);
    printf_solution();
    search();//BB搜索
    printf("read time: %f; init time: %f; search time: %f\n", read_time, init_time, get_utime() - read_time - init_time);
    printf("total time: %f\n", get_utime());
    //cout << ct << endl;
   // cout << LB << endl;
    //cout << nn << endl;
    printf("------------------------------------------------------------------------\n");
    freeAll();
    return 0;
}
