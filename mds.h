#ifndef _MDS_H
#include "Array.h"

struct Edge//边定义
{
    int v1,v2;
};

enum State//顶点状态定义
{
    Candidate = 0,          //剩余顶点
    Fixed,                            //固定，在最优解中
    Delete,                         //可删除节点
    Forbid,                         //禁止加入最优解，但又不可删除
};

struct Vertex//顶点数据结构
{
    State state;                //顶点状态
    int num_in_c;           //被支配的次数
    int cost;                       //顶点花费
    int score;						//加入该顶点后，自己及邻居从未支配到支配的顶点数量
    int is_in_search;       //是否在搜索集
    int degree;
    int *neighbor;
};

struct Vertex_sort//在顶点排序的时候方便记录顶点编号（用于顶点排序）
{
    int index;
    int degree;
};

Vertex *vertex;//顶点数组

class NeighborSet//闭邻居集合类
{
private:
    int _v;//顶点
    Array* neighbors;//邻居顶点集合
    int _validNeighborCnt = 0;//有效邻居个数，即未支配邻居个数
public:
    NeighborSet(int v, int capacity);
    ~NeighborSet();
    bool is_in_set(int);//判断一个顶点是否在闭邻居集合
    void add_neighbor(int);//添加邻居顶点
    int get_valid_neighbor_cnt();//获取未支配邻居数量
    int get_neighbor_cnt();//获取邻居总数
};

NeighborSet::NeighborSet(int v, int capacity)
{
     _v = v;
     neighbors = new Array(capacity);
}

NeighborSet::~NeighborSet()
{
    delete neighbors;
    neighbors = nullptr;
}

bool NeighborSet::is_in_set(int _v)
{
    return neighbors->is_in_array(_v);
}

void NeighborSet::add_neighbor(int _v)
{
    neighbors->insert_element(_v);
    if (vertex[_v].num_in_c == 0)
    {
        _validNeighborCnt++;
    }
}

int NeighborSet::get_valid_neighbor_cnt()
{
    return _validNeighborCnt;
}

int NeighborSet::get_neighbor_cnt()
{
    return neighbors->size();
}

#endif // !_MDS_H
