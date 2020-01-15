#ifndef _MDS_H
#include "Array.h"
#include <vector>

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
    int degree;                 //顶点度
    int *neighbor;          //邻居
};

struct Vertex_sort//在顶点排序的时候方便记录顶点编号（用于顶点排序）
{
    int index;
    int degree;
};

struct vertex_sort_score//用于顶点排序根据score值排序
{
    int index;
    int score;
};

// struct vertex_sort_valid_neighbor//用于子集缩减中有效闭邻居大小的比较
// {
//     int index;
//     int valid_neighbor_cnt;
//     std::vector<int> neighbors;
// };

Vertex *vertex;//顶点数组

class NeighborSet//闭邻居集合类
{
private:
    int _v;//顶点
    // int* _neighbors;//邻居顶点集合
    std::vector<int> _neighbors;
    int _validNeighborCnt = 0;//有效邻居个数，即未支配邻居个数
    int _size = 0;//容量
    int _capacity = 0;//容量上限
public:
    NeighborSet(int v, int capacity);
    ~NeighborSet();
    bool is_in_set(int);//判断一个顶点是否在闭邻居集合
    void add_neighbor(int);//添加邻居顶点
    int get_neighbor(int);//获取指定位置邻居
    int get_valid_neighbor_cnt();//获取未支配邻居数量
    int get_neighbor_cnt();//获取邻居总数
    int getV();//获取顶点V
};

NeighborSet::NeighborSet(int v, int capacity)
{
     _v = v;
     _capacity = capacity;
    //  _neighbors = new int[capacity];
    _neighbors = std::vector<int>(capacity);
}

NeighborSet::~NeighborSet()
{
    // delete[] _neighbors;
    // _neighbors = nullptr;
}

bool NeighborSet::is_in_set(int v)
{
    for (size_t i = 0; i < _size; i++)
    {
        if (_neighbors[i] == v)
            return true;
    }
    return false;
}

void NeighborSet::add_neighbor(int v)
{
    _neighbors[_size++] = v;
    if (vertex[v].num_in_c == 0)
    {
        _validNeighborCnt++;
    }
}

int NeighborSet::get_neighbor(int pos)
{
    return _neighbors[pos];
}

int NeighborSet::get_valid_neighbor_cnt()
{
    return _validNeighborCnt;
}

int NeighborSet::get_neighbor_cnt()
{
    return _size;
}

int NeighborSet::getV()
{
    return _v;
}

//从二维数组获取顶点，即邻居顶点
int getVertex(int i, int j)
{
    return vertex[i].neighbor[j];
}

#endif // !_MDS_H
