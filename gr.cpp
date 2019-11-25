#include <iomanip>
#include "gr.h"

long current_sol;

void buger()
{
    printf(" *\n");
}

void init() {
    std::cout << "initializing.." << endl;
    best_value = LONG_MAX; //取最大值
    cs_size = 0;
    current_sol = 0;
    t_length = 0;
    uncover_num = 0;
    //memset(t_index, -1, sizeof(t_index));
    for (size_t i = 0; i < vertex_num; i++)
    {
        t_index[i] = -1;
        uncover_vertex_index[i] = i;
        uncover_vertex[i] = i;
        uncover_num++;
    }
}

void simple_lock(int v)
{
    if (cs[v].locked == 1)
        return;
    lock_vertex(v, 1);
}

void init_reduce()
{
    int i, j;
    int a, b, c, sum;
    int v_neighbor, neighbor_num;
    memset(reduce, 0, sizeof(reduce));
    times(&start);
    std::cout << "first step--->: init reduce" << endl;
    for (i = 0; i < vertex_num; i++)
    {
        if (vertex_neightbourNum[i] == 1)
        {
            v_neighbor = vertex[i][0];
            if (cs[i].locked == 0)
            {
                reduce[i] = 1;
                //cs[v_neighbor].locked = 1;
                simple_lock(v_neighbor);
            }
        }
        else if (vertex_neightbourNum[i] == 2)
        {
            if (vertex_neightbourNum[vertex[i][0]] == 2)
            {
                a = vertex[i][0];
                b = vertex[i][1];
            }
            else if (vertex_neightbourNum[vertex[i][1]] == 2)
            {
                a = vertex[i][1];
                b = vertex[i][0];
            }
            else
                continue;
            if (vertex[a][0] == i)
                c = vertex[a][1];
            else
                c = vertex[a][0];
            if (b == c)
            {
                if (cs[i].locked == 0 && cs[a].locked == 0)
                    //if(cs[i].cost > cs[b].cost && cs[a].cost > cs[b].cost)
                {
                    reduce[i] = 1;
                    reduce[a] = 1;
                    simple_lock(b);
                }
            }
        }
        else if (vertex_neightbourNum[i] == 0)
        {
            simple_lock(i);
        }
        else
        {
            neighbor_num = vertex_neightbourNum[i];
            sum = 0;
            for (j = 0; j < neighbor_num; j++)
                if (vertex_neightbourNum[vertex[i][j]] == 1)
                    sum += cs[vertex[i][j]].cost;
            if (sum > cs[i].cost)
            {
                for (j = 0; j < neighbor_num; j++)
                {
                    if (vertex_neightbourNum[vertex[i][j]] == 1)
                    {
                        reduce[vertex[i][j]] = 1;
                    }
                }
                simple_lock(i);
            }
        }
    }
    int kkk, k;
    int lenn;
    remain_num = 0;
    for (i = 0; i < vertex_num; i++)
    {
        if (reduce[i] == 1)
        {
            //从图中删除该节点
            neighbor_num = vertex_neightbourNum[i];
            for (j = 0; j < neighbor_num; j++)
            {
                v_neighbor = vertex[i][j];
                lenn = vertex_neightbourNum[v_neighbor];
                for (k = 0; k < lenn; k++)
                    if (vertex[v_neighbor][k] == i)
                        break;
                vertex_neightbourNum[v_neighbor]--;
                kkk = vertex[v_neighbor][vertex_neightbourNum[v_neighbor]];
                vertex[v_neighbor][vertex_neightbourNum[v_neighbor]] = vertex[v_neighbor][k];
                vertex[v_neighbor][k] = kkk;
            }
        }
        else
        {
            //uncover_vertex[remain_num] = i;
            //uncover_vertex_index[i] = remain_num;
            remain_vertex[remain_num++] = i;
        }
    }
    uncover_num = remain_num;
    print_reduce_graph();
    superset_reduce();
}

//顶点从未支配中移除
void remove_from_uncover(int v)
{
    uncover_num--;
    int uk = uncover_vertex_index[v];
    int ck = uncover_vertex[uncover_num];
    uncover_vertex[uk] = ck;
    uncover_vertex_index[ck] = uk;
}

//顶点score值改变，判断加入t还是移出t
void modify_t(int v)
{
    if (cs[v].score == 0)
        reduce[v] = 1;
    if (cs[v].score == 0 && cs[v].is_in_c == 1 && t_index[v] == -1 && cs[v].locked == 0)
    {
        t[t_length] = v;
        t_index[v] = t_length++;
    }
    else if (t_index[v] != -1 && cs[v].score != 0)
    {
        t_length--;
        int uk = t_index[v];
        int ck = t[t_length];
        t[uk] = ck;
        t_index[ck] = uk;
        t_index[v] = -1;
    }
}

//把顶点c固定
void lock_vertex(int c, int locked_add)
{
    cs[c].locked = 1;
    cs[c].is_in_c = 1;
    if (locked_add == 1)
    {
        cs_vertex[cs_size] = c;
        cs_vertex_index[c] = cs_size++;
    }
    cs[c].score = -cs[c].score;
    if (cs[c].score == 0 && t_index[c] == -1 && cs[c].locked == 0)
    {
        t[t_length] = c;
        t_index[c] = t_length++;
    }
    //if (cs[c].num_in_c == 0)
    //	remove_from_uncover(c);
    for (size_t h = 0; h < vertex_neightbourNum[c]; h++)
    {
        int v_n = vertex[c][h];
        //if (cs[i].num_in_c == 0)
        //	remove_from_uncover(i);
        cs[v_n].num_in_c++;
        if (cs[c].num_in_c == 0)
        {
            //v_n被支配，修改score
            //c 支配次数 0->1
            cs[v_n].score -= vertex_weight[c];
            modify_t(v_n);
        }
        else if (cs[c].num_in_c == 1 && cs[v_n].is_in_c == 1)
        {
            //v_n存在在候选解, 此时v_n的score值为负
            //c 支配次数 1->2
            cs[v_n].score += vertex_weight[c];
            modify_t(v_n);
        }
        if (cs[v_n].is_in_c == 1)
        {
            if (cs[v_n].num_in_c == 2)
                //v_n已经在候选解，c加入候选解后，v_n又被支配了一次
                cs[v_n].score += vertex_weight[v_n];
            modify_t(v_n);
            continue;
        }
        //处理二层邻居
        int cnt = 0;
        int s = 0;
        for (size_t l = 0; l < vertex_neightbourNum[v_n]; l++)
        {
            int j = vertex[v_n][l];
            if (j == c)
                continue;
            if (cs[j].is_in_c)
            {
                s = j;
                cnt++;
            }
        }
        if (cs[v_n].is_in_c)
        {
            s = v_n;
            cnt++;
        }
        if (cnt == 0)
        {
            //c是i的邻居里面第一个加入候选解的
            cs[v_n].score -= vertex_weight[v_n];
            for (size_t l = 0; l < vertex_neightbourNum[v_n]; l++)
            {
                int j = vertex[v_n][l];
                if (j == c)
                    continue;
                //i被c支配，所以二层邻居j score值减少
                cs[j].score -= vertex_weight[v_n];
                modify_t(j);
            }
        }
        else if (cnt == 1)
        {
            //c是这里组里第二个加入候选解的顶点
            cs[s].score += vertex_weight[v_n];
            modify_t(s);
        }
        modify_t(v_n);
    }
    cs[c].num_in_c++;
}

void superset_reduce()
{
    int iter = 0;
    std::cout << "second step--->: superset reduce" << endl;
    queue<int> q_searchset;
    for (size_t i = 0; i < vertex_num; i++)
        if (cs[i].num_in_c == 0)
        {
            cs[i].is_in_search = 1;
            q_searchset.push(i);
        }
    while (!q_searchset.empty()) {
        //时间计算
        times(&finish);
        double tt = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK);
        tt = round(tt * 100) / 100.0;
        if (time_limit != 0 && tt > time_limit)
            break;
        iter++;
        int v = q_searchset.front();
        q_searchset.pop();
        cs[v].is_in_search = 0;
        if (cs[v].num_in_c != 0)
            continue;
        int set_count = 0;
        //存放v闭邻居score>0集合
        vector<int> node;
        vector<int> node_score;
        node.push_back(v);
        node_score.push_back(cs[v].score);
        set_count++;
        for (int j = 0; j < vertex_neightbourNum[v]; ++j) {
            int v_neighbor = vertex[v][j];
            if (cs[v_neighbor].score > 0) {
                node.push_back(v_neighbor);
                node_score.push_back(cs[v_neighbor].score);
                ++set_count;
            }
        }
        if (set_count == 1)
            lock_vertex(v, 1);
        else {
            //寻找score最大的节点
            int max_score_index = 0;
            int max_score = node_score[0];
            for (int i = 1; i < set_count; ++i) {
                if (node_score[i] > max_score) {
                    max_score = node_score[i];
                    max_score_index = i;
                }
            }
            int max_score_v = node[max_score_index];
            auto max_set = new NeighborSet(max_score_v);
            int cnt = 0;
            //初始化超集候选集合元素
            if (cs[max_score_v].num_in_c == 0)
                max_set->neighbors[cnt++] = max_score_v;
            for (int i = 0; i < vertex_neightbourNum[max_score_v] && cnt < max_score; ++i) {
                int v_neighbor = vertex[max_score_v][i];
                if (cs[v_neighbor].num_in_c == 0)
                    max_set->neighbors[cnt++] = v_neighbor;
            }
            //判断score值最大的集合是不是其他集合的超集
            bool is_super_set = true;
            for (int i = 0; i < set_count; ++i) {
                if (node[i] == max_score_v)
                    continue;
                int v_neighbor = node[i];
                int v_neighbor_score = node_score[i];
                cnt = 0;
                if (cs[v_neighbor].num_in_c == 0)
                {
                    if (!max_set->is_in_set(v_neighbor))
                    {
                        is_super_set = false;
                        break;
                    }
                    else
                        cnt++;
                }
                for (int j = 0; j < vertex_neightbourNum[v_neighbor] && cnt < v_neighbor_score; ++j) {
                    int v_neighbor_neighbor = vertex[v_neighbor][j];
                    if (cs[v_neighbor_neighbor].num_in_c == 0) {
                        if (max_set->is_in_set(v_neighbor_neighbor)) {
                            cnt++;
                        } else {
                            is_super_set = false; //不是超集，跳出循环
                            break;
                        }
                    }
                }
                if (!is_super_set)
                    break;
            }
            if (is_super_set)
            {
                lock_vertex(max_set->v, 1);
                //max_set->v的所有未覆盖的二层邻居都加入队列
                for (size_t i = 0; i < vertex_neightbourNum[max_set->v]; i++) {
                    int v_n = vertex[max_set->v][i];
                    for (size_t j = 0; j < vertex_neightbourNum[v_n] ; j++) {
                        int v_n_n = vertex[v_n][j];
                        //v_n_n邻居有未被支配的顶点
                        //v_n的支配次数不止一次，也就是在加入顶点max_set->v后，v_n_n的score值发生变化
                        if (cs[v_n_n].score > 0 && cs[v_n].num_in_c == 1) 
                        {
                            int cnt2 = 0;
                            if (cs[v_n_n].num_in_c == 0 && cs[v_n_n].is_in_search == 0)
                            {
                                cs[v_n_n].is_in_search = 1;
                                q_searchset.push(v_n_n);
                                cnt2++;
                            }
                            for (size_t k = 0; k < vertex_neightbourNum[v_n_n] && cnt2 < cs[v_n_n].score; k++)
                            {
                                int v_n_n_n = vertex[v_n_n][k];
                                if (cs[v_n_n_n].num_in_c == 0 && cs[v_n_n_n].is_in_search == 0)
                                {
                                    cs[v_n_n_n].is_in_search = 1;
                                    q_searchset.push(v_n_n_n);
                                    cnt2++;
                                }
                            }
                        }
                    }
                }
            }
            delete max_set;
        }
    }
    cout << "iter: " << iter << endl;
    // print_reduce_graph();
    // print_density();
    // print_degree();
    // times(&finish);
    // double tt = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
    // tt= round(tt * 100)/100.0;
    // cout << "Time: " << tt << "s" <<endl;
    subset_reduce();
}

void subset_reduce()
{
    cout << "third step--->:subset reduce" <<endl;
    for (size_t i = 0; i < vertex_num; i++)
    {
        if (cs[i].num_in_c == 0)
        {
            int cnt = 1;    //i本身也要存在集合里面
            for (size_t j = 0; j < vertex_neightbourNum[i]; j++)
            {
                int v_n = vertex[i][j];
                if (cs[v_n].score > 0)
                {
                    cnt++;
                }
            }
            //构造所有邻居集合
            vector<NeighborSet*> sets(cnt);
            sets[0] = new NeighborSet(i, cnt);
            sets[0]->addNeighbor(i);
            int cnt1 = 0;
            for (size_t j = 0; j < vertex_neightbourNum[i] && cnt1 < cnt; j++)
            {
                int v_n = vertex[i][j];
                if (cs[v_n].score > 0)
                {
                    sets[0]->addNeighbor(v_n);
                    cnt1++;
                    sets[cnt1] = new NeighborSet(v_n);
                    int cnt2 = 0;
                    if (cs[v_n].num_in_c == 0)
                    {
                        sets[cnt1]->addNeighbor(v_n);
                        cnt2++;
                    }
                    for (size_t k = 0; k < vertex_neightbourNum[v_n] && cnt2 < cs[v_n].score; k++)
                    {
                        int v_n_n = vertex[v_n][k];
                        if (cs[v_n_n].num_in_c == 0)
                        {
                            sets[cnt1]->addNeighbor(v_n_n);
                            cnt2++;
                        }
                    }
                }
            }
            //两层循环遍历删除子集顶点
            for (size_t j = 0; j < sets.size(); j++)
            {
                for (size_t k = j + 1; k < sets.size(); k++)
                {
                    auto set_a = sets[j];
                    auto set_b = sets[k];
                    //外层节点已经被删除，打断内部for循环
                    if (reduce[set_a->v] == 1)
                    {
                        break;
                    }
                    //内部节点已经被删除，跳过本次内部循环
                    if (reduce[set_b->v] == 1)
                    {
                        continue;
                    }
                    //如果set_b元素数量比set_a多，交换指针，set_b始终是待判断子集的指针
                    if (set_a->neighbor_cnt < set_b->neighbor_cnt)
                    {
                        auto temp = set_b;
                        set_b = set_a;
                        set_a = temp;
                    }
                    //判断set_b是不是set_a的子集
                    bool is_subset = true;
                    for (size_t t = 0; t < set_b->neighbor_cnt; t++)
                    {
                        if (!set_a->is_in_set(set_b->neighbors[t]))
                        {
                            is_subset = false;
                            break;
                        }
                    }
                    if (is_subset)
                    {
                        //set_b是set_a的子集，把set_b的头元素即v删除
                        reduce[set_b->v] = 1;
                    }
                }
            }
            //free memory
            for (size_t j = 0; j < sets.size(); j++)
            {
                delete sets[j];
            }
        }
    }
    print_reduce_graph();
    print_density();
    print_degree();
    times(&finish);
    double tt = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
    tt= round(tt * 100)/100.0;
    cout << "Time: " << tt << "s" <<endl;
}

void print_reduce_graph()
{
    int locked_num = 0;
    int uncover_num = 0;
    int remove_num = 0;
    remain_num = 0;
    int score_lock_0 = 0;
    for (size_t i = 0; i < vertex_num; i++)
    {
        if  (cs[i].locked == 1)
        {
            locked_num++;
            if (cs[i].score == 0)
            {
                score_lock_0++;
            }
            continue;
        }
        if (reduce[i] == 1)
        {
            remove_num++;
            continue;
        }
        remain_num++;
        if (cs[i].num_in_c == 0)
        {
            uncover_num++;
        }
    }
//    cout << "score_lock_0: " << score_lock_0 << endl;
//    if (remain_num + locked_num + remove_num == vertex_num)
//        cout << "check ok" << endl;
//    else
//        cout << vertex_num - remain_num - locked_num - remove_num << endl;
    std::cout << "Total Vertex: " << vertex_num << endl;
    std::cout << "Delete Vertex: " << remove_num << endl;
    std::cout << "Fixed Vertex: " << locked_num << endl;
    std::cout << "Uncover Vertex: " << uncover_num << endl;
    std::cout << "Remain Vertex:" << remain_num << endl;
    std::cout << "Percent: " <<  fixed << setprecision(2) << remain_num * 1.0 / vertex_num * 100 << "%" << endl;
}

void print_density()
{
    int edge_cnt = 0;
    std::ofstream openfile("abc", std::ios::out);
    for (size_t i = 0; i < vertex_num; i++)
    {
        if (cs[i].locked == 1 || reduce[i] == 1)
            continue;
        for (int j = 0; j < vertex_neightbourNum[i]; ++j) {
            int v_n = vertex[i][j];
            if (v_n > i && cs[v_n].locked != 1 && reduce[v_n] != 1)
            {
                edge_cnt++;
                openfile << "e " << i << " " << v_n <<endl;
            }
        }

//        if (cs[i].score == 0 || cs[i].locked == 1)
//        {
//            continue;
//        }
//        for (size_t j = 0; j < vertex_neightbourNum[i]; j++)
//        {
//            int v_n = vertex[i][j];
//            if (v_n > i && cs[v_n].score != 0 && cs[v_n].locked != 1)
//            {
//                edge_cnt++;
//            }
//        }
    }
    openfile.close();
    cout << "Remain vertex: " << remain_num << endl;
    cout << "Remain edges: " << edge_cnt << endl;
    cout << "Density1: " << fixed << setprecision(3) << edge_num *1.0 / vertex_num / (vertex_num - 1) * 100 << "%" << endl;
    if (remain_num != 0)
    {
        cout << "Density2: " << fixed << setprecision(3) << edge_cnt * 1.0 / remain_num / (remain_num - 1) * 100 << "%" << endl;
    }
    else
    {
        cout << "Density2: 0" << endl;
    }
}

void print_degree()
{
    int d1 = 0, d2 = 0, d3 = 0;
    for (size_t i = 0; i < vertex_num; i++)
    {
        if (cs[i].score == 0 || cs[i].locked == 1)
        {
            continue;
        }
        int n_cnt = 0;
        for (size_t j = 0; j < vertex_neightbourNum[i]; j++)
        {
            int v_n = vertex[i][j];
            if (cs[v_n].score != 0 && cs[v_n].locked != 1)
            {
                n_cnt++;
            }
        }
        if (n_cnt == 1)
        {
            d1++;
        }
        else if (n_cnt == 2)
            d2++;
        else if (n_cnt == 3)
            d3++;
    }
    cout << "d1:" << d1 << endl;
    cout << "d2:" << d2 << endl;
    cout << "d3:" << d3 << endl;
}

inline int compare(int s1, int c1, int s2, int c2){
    if(c1==c2) {
        if(s1>s2) return 1;
        else if(s1==s2) return 0;
        else return -1;
    }
    long long t1=s1, t2=s2;
    t1=t1*c2;
    t2=t2*c1;
    if(t1>t2) return 1;
    else if(t1==t2) return 0;
    else return -1;
}

int check(){ // check if the solution is a correct cover
    int i,j,k,v;
    for(v=0;v<vertex_num;v++)
    {
        reduce[v] = 0;
    }
    for(v=0;v<vertex_num;v++)
    {
        j = v;
        if(best_sol[j]==0 && cs[j].locked == 1)
            return 0;
        if(best_sol[j]==0)
            continue;
        reduce[j]=1;
        for(k=0;k<vertex_neightbourNum_bak[j];k++)
        {
            reduce[vertex[j][k]]=1;
        }
    }
    for(v=0;v<vertex_num;v++){
        i = v;
        if(!reduce[i])
            return 0;
    }
    return 1;
}

void add(int c, int locked_add, int init_add){
    int i,j,k,cnt,s, ii, jj, ix,h, l;
    int uk,ck;
    cs[c].is_in_c=1;
    current_sol += cs[c].cost;
    if(locked_add == 1)
    {
        cs_vertex[cs_size] = c;
        cs_vertex_index[c] = cs_size;
        cs_size++;
    }

    cs[c].score=-cs[c].score;
    if(cs[c].score == 0 && t_index[c] == -1 && cs[c].locked == 0)
    {
        t[t_length] = c;
        t_index[c] = t_length;
        t_length++;
    }

    if( cs[c].num_in_c==0)
    {
        uncover_num--;
        if(init_add == 0)
        {
            uk = uncover_vertex_index[c];
            ck = uncover_vertex[uncover_num];
            uncover_vertex[uk] = ck;
            uncover_vertex_index[ck] = uk;
        }
    }
    for(h=0;h<vertex_neightbourNum[c]; h++)
    {//C集合中每一个变量h   处理未覆盖集合
        i=vertex[c][h];

        if(cs[i].num_in_c==0)
        {
            uncover_num--;
            if(init_add == 0)
            {
                uk = uncover_vertex_index[i];
                ck = uncover_vertex[uncover_num];
                uncover_vertex[uk] = ck;
                uncover_vertex_index[ck] = uk;
            }
        }

        cs[i].num_in_c++;/////////////////////
        cnt=0;
        if( cs[c].num_in_c==0){
            cs[i].score-=vertex_weight[c];
        }    else if( cs[c].num_in_c==1&& cs[i].is_in_c==1){
            cs[i].score+=vertex_weight[c];
        }


        if(cs[i].is_in_c && init_add == 0)
        {
            if(cs[i].num_in_c == 2)
                cs[i].score += vertex_weight[i];
            if(cs[i].score == 0 && t_index[i] == -1 && cs[i].locked == 0)
            {
                t[t_length] = i;
                t_index[i] = t_length;
                t_length++;
            }
            else if(t_index[i] != -1 && cs[i].score != 0)
            {
                t_length--;
                uk = t_index[i];
                ck = t[t_length];
                t[uk] = ck;
                t_index[ck] = uk;
                t_index[i] = -1;
            }
            continue;
        }


        for(l=0;l<vertex_neightbourNum[i];l++)
        {//因为这个变量，变成这个集合的邻居集合
            j=vertex[i][l];

            if(j==c) continue;
            if(cs[j].is_in_c)
            {
                s=j;
                cnt++;
            }
        }
        if(cs[i].is_in_c)
        {
            s=i;
            cnt++;
        }

        if(cnt==0){ // c is the first one covering this row in C
            cs[i].score-=vertex_weight[i];////////////////
            for(l=0; l<vertex_neightbourNum[i]; l++){
                j=vertex[i][l];
                if(j==c)
                    continue;
                cs[j].score-=vertex_weight[i];//候选解中不覆盖这个变量，所以当覆盖之后，所以以前覆盖这个变量的集合score取值必须减去这个变量
                if(cs[j].score == 0 && cs[j].is_in_c == 1 && t_index[j] == -1 && cs[j].locked == 0)
                {
                    t[t_length] = j;
                    t_index[j] = t_length;
                    t_length++;
                }
                else if(t_index[j] != -1 && cs[j].score != 0)
                {
                    t_length--;
                    uk = t_index[j];
                    ck = t[t_length];
                    t[uk] = ck;
                    t_index[ck] = uk;
                    t_index[j] = -1;
                }
            }
        } else if(cnt==1){// c is second one covering this row in C
            cs[s].score+=vertex_weight[i];//候选解中覆盖这个变量一次，所以加入这个集合以后，所以以前覆盖这个变量的集合score取值必须加上这个变量
            if(cs[s].score == 0 && cs[s].is_in_c == 1 && t_index[s] == -1 && cs[s].locked == 0)
            {
                t[t_length] = s;
                t_index[s] = t_length;
                t_length++;
            }
            else if(t_index[s] != -1 && cs[s].score != 0)
            {
                t_length--;
                uk = t_index[s];
                ck = t[t_length];
                t[uk] = ck;
                t_index[ck] = uk;
                t_index[s] = -1;
            }
        }

        if(cs[i].score == 0 && cs[i].is_in_c == 1 && t_index[i] == -1 && cs[i].locked == 0)
        {
            t[t_length] = i;
            t_index[i] = t_length;
            t_length++;
        }
        else if(t_index[i] != -1 && cs[i].score != 0)
        {
            t_length--;
            uk = t_index[i];
            ck = t[t_length];
            t[uk] = ck;
            t_index[ck] = uk;
            t_index[i] = -1;
        }
    }
    cs[c].num_in_c++;
}

void remove(int c, int init_remove){
    int uk, ck;
    cs[c].is_in_c=0;
    current_sol -= cs[c].cost;
    cs_size--;
    int kn = cs_vertex_index[c];
    int bn = cs_vertex[cs_size];
    cs_vertex[kn] = bn;
    cs_vertex_index[bn] = kn;

    cs[c].score=-cs[c].score;
    if(cs[c].score == 0 && t_index[c] != -1)
    {
        t_length--;
        uk = t_index[c];
        ck = t[t_length];
        t[uk] = ck;
        t_index[ck] = uk;
        t_index[c] = -1;
    }

    if( cs[c].num_in_c==1)
    {
        //uncover_vertex[uncover_num] = c;//////////////////////////
        // uncover_vertex_index[c] = uncover_num;
        uncover_num++;
    }

    int i,j,k,cnt,s,h,l;
    for(h=0; h<vertex_neightbourNum[c]; h++){
        i=vertex[c][h];
        if(cs[i].num_in_c==1)
        {
            //uncover_vertex[uncover_num] = i;
            //uncover_vertex_index[c] = uncover_num;
            //uncover_vertex_index[i] = uncover_num;
            uncover_num++;
        }
        cs[i].num_in_c--;
        cnt=0;
        if( cs[c].num_in_c==2&& cs[i].is_in_c==1){
            cs[i].score-=vertex_weight[c];
        }    else if( cs[c].num_in_c==1){
            cs[i].score+=vertex_weight[c];
        }


        if(cs[i].is_in_c && init_remove == 0)
        {
            if(cs[i].num_in_c == 1)
                cs[i].score -= vertex_weight[i];
            if(cs[i].score == 0 && t_index[i] == -1 && cs[i].locked == 0)
            {
                t[t_length] = i;
                t_index[i] = t_length;
                t_length++;
            }
            else if(t_index[i] != -1 && cs[i].score != 0)
            {
                t_length--;
                uk = t_index[i];
                ck = t[t_length];
                t[uk] = ck;
                t_index[ck] = uk;
                t_index[i] = -1;
            }
            continue;
        }
        for(l=0;l<vertex_neightbourNum[i];l++){
            j=vertex[i][l];
            if(j==c) continue;
            if(cs[j].is_in_c){
                cnt++;
                s=j;
            }
        }
        if(cs[i].is_in_c){
            s=i;
            cnt++;
        }
        if(cnt==0){
            cs[i].score+=vertex_weight[i];////////////////
            for(l=0;l<vertex_neightbourNum[i];l++){
                j=vertex[i][l];
                if(j==c)
                    continue;
                cs[j].score+=vertex_weight[i];
                if(cs[j].score == 0 && cs[j].is_in_c == 1 && t_index[j] == -1 && cs[j].locked == 0)
                {
                    t[t_length] = j;
                    t_index[j] = t_length;
                    t_length++;
                }
                else if(t_index[j] != -1 && cs[j].score != 0)
                {
                    t_length--;
                    uk = t_index[j];
                    ck = t[t_length];
                    t[uk] = ck;
                    t_index[ck] = uk;
                    t_index[j] = -1;
                }
            }
        } else if(cnt==1){
            cs[s].score-=vertex_weight[i];
            if(cs[s].score == 0 && cs[s].is_in_c == 1 && t_index[s] == -1 && cs[s].locked == 0)
            {
                t[t_length] = s;
                t_index[s] = t_length;
                t_length++;
            }
            else if(t_index[s] != -1 && cs[s].score != 0)
            {
                t_length--;
                uk = t_index[s];
                ck = t[t_length];
                t[uk] = ck;
                t_index[ck] = uk;
                t_index[s] = -1;
            }
        }
        if(cs[i].score == 0 && cs[i].is_in_c == 1 && t_index[i] == -1 && cs[i].locked == 0)
        {
            t[t_length] = i;
            t_index[i] = t_length;
            t_length++;
        }
        else if(t_index[i] != -1 && cs[i].score != 0)
        {
            t_length--;
            uk = t_index[i];
            ck = t[t_length];
            t[uk] = ck;
            t_index[ck] = uk;
            t_index[i] = -1;
        }
    }
    cs[c].num_in_c--;

}

void uncov_r_weight_inc(){
    int i,j,h,v,nn;
    for(v = 0; v < uncover_length; v++)
    {
        i = uncover_vertex[v];
        if(cs[i].num_in_c == 0)
        {
            vertex_weight[i] += 1;
            cs[i].score += 1;
            nn = vertex_neightbourNum[i];
            for(h = 0; h < nn; h++)
            {
                j = vertex[i][h];
                cs[j].score += 1;
            }
        }
    }
}

void update_best_sol(){
    int i,j,v;
    if(current_sol < best_value){
        best_value = current_sol;
        for(v=0;v<remain_num;v++){
            j = remain_vertex[v];
            best_sol[j]=0;
            if(cs[j].is_in_c){
                best_sol[j]=1;
            }
        }

        times(&finish);
        real_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_time = round(real_time * 100)/100.0;
    }
}

int main(int argc, char *argv[]){
    if(argc<2){
        printf("input wrong\n");
        return 0;
    }

    build_instance_massive(argv[1]);
    time_limit=atof(argv[2]);

    init();
    init_reduce();
    free_all();
    return 0;
}
