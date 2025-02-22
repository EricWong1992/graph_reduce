#include "cc42.h"

int bms_thre = 1;
int bms;
int each_iter=0;
long current_sol;
int remove_num1, remove_num2;


void buger()
{
    printf(" *\n");
}

void init(){
    best_value = LONG_MAX; //取最大值
    cs_size = 0;
    locked_num = 0;
    tabu_list.clear();
    current_sol = 0;
    t_length = 0;
   //memset(t_index, -1, sizeof(t_index));
    for(int i = 0; i < vertex_num; i++)
    {
        t_index[i] = -1;
    }
}

void init_reduce()
{
    int i,j,v,kn;
    int a,b,c,sum;
    int v_neighbor, neighbor_num;
    memset(reduce, 0, sizeof(reduce));
    for(i = 0; i < vertex_num; i++)
    {
        if(vertex_neightbourNum[i] == 1)
        {
            v_neighbor = vertex[i][0];
            //if(cs[i].cost > cs[v_neighbor].cost)
            //if(vertex_neightbourNum[v_neighbor] > 1)
            if(cs[i].locked == 0)
            {
                reduce[i] = 1;
                cs[v_neighbor].locked = 1;
            }
        }
        else if(vertex_neightbourNum[i] == 2)
        {
            if(vertex_neightbourNum[vertex[i][0]] == 2)
            {
                a = vertex[i][0];
                b = vertex[i][1];
            }
            else if(vertex_neightbourNum[vertex[i][1]] == 2)
            {
                a = vertex[i][1];
                b = vertex[i][0];
            }
            else
                continue;
            if(vertex[a][0] == i)
                c = vertex[a][1];
            else
                c = vertex[a][0];
            if(b == c)
            {
                if(cs[i].locked == 0 && cs[a].locked == 0)
                //if(cs[i].cost > cs[b].cost && cs[a].cost > cs[b].cost)
                {
                    reduce[i] = 1;
                    reduce[a] = 1;
                    cs[b].locked = 1;
                }
            }
        }
        else if(vertex_neightbourNum[i] == 0)
        {
            cs[i].locked = 1;
        }
        else
        {
            neighbor_num = vertex_neightbourNum[i];
            sum = 0;
            for(j = 0; j < neighbor_num; j++)
                if(vertex_neightbourNum[vertex[i][j]] == 1)
                    sum += cs[vertex[i][j]].cost;
            if(sum > cs[i].cost)
            {
                for(j = 0; j < neighbor_num; j++)
                    if(vertex_neightbourNum[vertex[i][j]] == 1)
                    {
                        reduce[vertex[i][j]] = 1;
                    }
                cs[i].locked = 1;
            }
        }
    }
    int kkk,k;
    int lenn;
    remain_num = 0;
    for(i = 0; i < vertex_num; i++)
    {
            if(reduce[i])
            {
                neighbor_num = vertex_neightbourNum[i];
                for(j = 0; j < neighbor_num; j++)
                {
                    v_neighbor = vertex[i][j];
                    lenn = vertex_neightbourNum[v_neighbor];
                    for(k = 0; k < lenn; k++)
                            if(vertex[v_neighbor][k] == i)
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
    //uncover_num = remain_num;
   // printf("%d %d ", vertex_num, remain_num);
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
/*
int check(){ // check if the solution is a correct cover
    int i,j,k,v;
    for(v=0;v<remain_num;v++)
    {
            j = remain_vertex[v];
            reduce[j] = 0;
    }

    for(v=0;v<remain_num;v++){
            j = remain_vertex[v];
            if(best_sol[j]==0 && cs[j].locked == 1)
                return 0;
     if(best_sol[j]==0) continue;
     reduce[j]=1;
     for(k=0;k<vertex_neightbourNum[j];k++){
       reduce[vertex[j][k]]=1;
     }
    }
    for(v=0;v<remain_num;v++){
            i = remain_vertex[v];
     if(!reduce[i]) return 0;
    }
    return 1;
}*/

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
        for(k=0;k<vertex_neightbourNum1[j];k++)
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


void init_best(){
    int cnt, kn;
    int i,j,k,l,jj,v;
    int sr, ct;
    for(v = 0; v < remain_num; v++)
    {
        i = remain_vertex[v];
        if(cs[i].locked)
        {
            cs[i].is_in_c = 1;
            cs[i].num_in_c++;
            current_sol += cs[i].cost;
            kn = vertex_neightbourNum[i];
            for(j = 0; j < kn; j++)
                cs[vertex[i][j]].num_in_c++;
        }
    }
    uncover_num = 0;
    for(v = 0; v < remain_num; v++)
    {
        i = remain_vertex[v];
        if(cs[i].locked)
            continue;
        kn = vertex_neightbourNum[i];
        for(k = 0; k <kn; k++)
        {
            j = vertex[i][k];
            if(cs[j].num_in_c == 0)
                cs[i].score += 1;
        }
        if(cs[i].num_in_c == 0)
        {
            cs[i].score += 1;
            uncover_vertex[uncover_num] = i;
            uncover_vertex_index[i] = uncover_num;
            uncover_num++;
        }
    }
   //
   //wcnf
   //
   //printf("%d  ", vertex_num - uncover_num);
   //printf("\n");
    free_all();
    exit(0);
    int sst = 0;
    int uncover_v;
    while(uncover_num>0)
    {
        /*
        if(sst % 2000 == 0)
        {
            times(&finish);
                        double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
                        finish_time = round(finish_time * 100)/100.0;
                printf("%d %d %f", sst , uncover_num, finish_time);
                buger();
        }
        sst++;*/
        uncover_v = uncover_vertex[rand()%uncover_num];

        sr = cs[uncover_v].score;
        ct = cs[uncover_v].cost;
        best_array[0] = uncover_v;
        cnt = 1;
        kn = vertex_neightbourNum[uncover_v];
        for(v=0;v<kn;v++){
            j = vertex[uncover_v][v];
            if(cs[j].is_in_c) continue;
            //if(cs[j].num_in_c > 0) continue;
            k=compare(sr,ct,cs[j].score, cs[j].cost);
            if(sr==INT_MIN||k<0){
                sr=cs[j].score;
                ct=cs[j].cost;
                best_array[0]=j;
                cnt=1;
            } 
            else if(k==0){
	            best_array[cnt++]=j;
            }
        }
        if(cnt>0){
            l=rand()%cnt;
            add(best_array[l], 1, 0);
        }
    }
    update_best_sol();
    if(check()==0){
        printf("initial wrong\n");exit(0);
    }
    printf("%ld	%.2f	", best_value, real_time);
}


void add(int c, int locked_add, int init_add){ 
    int i,j,k,cnt,s, ii, jj, ix;
    int uk,ck;
    cs[c].is_in_c=1;
    current_sol += cs[c].cost;
    if(locked_add == 1)
    {
        cs_vertex[cs_size] = c;
        cs_vertex_index[c] = cs_size++;
    }

    //刚添加到c的节点score值转负，score越小节点价值越大
    cs[c].score=-cs[c].score;
    if(cs[c].score == 0 && t_index[c] == -1 && cs[c].locked == 0)
    {
        //记录score为0的节点
        t[t_length] = c;
        t_index[c] = t_length++;
    }

    if (cs[c].num_in_c==0)
    {
        uncover_num--;
        if(init_add == 0)
        {
            //更新uncover集合
            uk = uncover_vertex_index[c];
            ck = uncover_vertex[uncover_num];
            uncover_vertex[uk] = ck;
            uncover_vertex_index[ck] = uk;
        }
    }
    for(size_t h=0;h<vertex_neightbourNum[c]; h++)
    {
        //C集合中每一个变量h   处理未覆盖集合
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

        cs[i].num_in_c++;   //未被支配的邻居现在已经支配
        if(cs[c].num_in_c==0)
        {
            //c没有邻居没有添加到候选解，对邻居i减去c weight
            //c 支配次数 0->1
            cs[i].score-=vertex_weight[c];
        }
        else if (cs[c].num_in_c==1 && cs[i].is_in_c==1)
        {
            //i已经在候选解，且i支配c，此时i的score是负
            //c 支配次数 1->2
            cs[i].score+=vertex_weight[c];
        }

        if(cs[i].is_in_c && init_add == 0)
        {
            if(cs[i].num_in_c == 2)
                //TODO？？？
                cs[i].score += vertex_weight[i];
            if(cs[i].score == 0 && t_index[i] == -1 && cs[i].locked == 0)
            {
                t[t_length] = i;
                t_index[i] = t_length++;
            }
            else if(t_index[i] != -1 && cs[i].score != 0)
            {
                uk = t_index[i];
                ck = t[--t_length];
                t[uk] = ck;
                t_index[ck] = uk;
                t_index[i] = -1;
            }
            continue;
        }

        cs[i].config = 1;
        //if(init_add == 0)
                //cs[i].config = 2;
        //二层邻居处理
        cnt = 0;
        for(size_t l=0;l<vertex_neightbourNum[i];l++)
        {//因为这个变量，变成这个集合的邻居集合
            j=vertex[i][l];

            if(j==c) continue;
            if(cs[j].is_in_c)
            {
                s=j;
                cnt++;
            }
            //if(rand()%100 < 50)
            //if(init_add == 1)
            cs[j].config = 2;
        }
        if(cs[i].is_in_c)
        {
            s=i;
            cnt++;
        }

        if(cnt==0)
        { 
            // c is the first one covering this row in C
            cs[i].score-=vertex_weight[i];////////////////
            for(size_t l=0; l<vertex_neightbourNum[i]; l++)
            {
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
        }
        else if(cnt==1)
        {// c is second one covering this row in C
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
    //if(init_remove == 1)
                cs[c].config=0;

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
        //if(init_remove == 1)
     cs[i].config = 2;
        for(l=0;l<vertex_neightbourNum[i];l++){
            j=vertex[i][l];
            if(j==c) continue;
            if(cs[j].is_in_c){
	cnt++;
	s=j;
            }
                //if(rand()%100 < 50)
                //if(init_remove == 1)
         cs[j].config=2;
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

int in_tabu(int i){
    return tabu_list.find(i)!=tabu_list.end();
}
int state=0;

int find_best_in_c_simp(int allowTabu){
    int i, maxc,j,k,v;

    int sr=INT_MIN, ct=1;
    int kn = cs_size;
    state = 0;
    for(v=0;v<kn;v++){
                i = v;
                i = cs_vertex[i];
                if(allowTabu&&in_tabu(i))
                        continue;
                state = 1;
                k=compare(sr,ct, cs[i].score, cs[i].cost);
                if(sr==INT_MIN||k<0){
                        sr=cs[i].score;
                        ct=cs[i].cost;
                        //maxc=i;
                        maxc=i;
                } else if(k==0){
                if(cs[maxc].time_stamp>cs[i].time_stamp){
                        //maxc=i;
                        maxc=i;
                }
                }
    }
    return maxc;
}

int find_best_in_c(int allowTabu){
    int i, maxc,j,k,v;
    int sr=INT_MIN, ct=1;
    int temp_int;
    int kn,flag;
    state=0;
    double r_n = rand()/(RAND_MAX+1.0);
    //double p = 1.0 / (double)best_sol_found;
    double p = exp(-each_iter);
    if(r_n < p)
    {
                //kn = cs_size;
                //bms = INT_MAX;
                kn = 1024;
    }
    else
    {
                /*bms = 100;
                temp_int = bms_thre % cs_size;
                if(bms < temp_int)
                        bms = temp_int;
                kn = cs_size;
                if(kn > bms)
                                kn = bms;*/
                kn = 50 + rand()%10;
    }
                if(kn > cs_size)
                        kn = cs_size;
                for(v=0;v<kn;v++){
                        if(cs_size == kn)
                                i = v;
                        else
                                i = rand()%cs_size;
                        i = cs_vertex[i];
                        if(allowTabu&&in_tabu(i)) continue;
                        state=1;
                        k=compare(sr,ct, cs[i].score, cs[i].cost);
                        if(sr==INT_MIN||k<0){
                                sr=cs[i].score;
                                ct=cs[i].cost;
                                maxc=i;
                        } else if(k==0){
                                        if(cs[maxc].time_stamp>cs[i].time_stamp){
                                                maxc=i;
                                        }
                        }
                }
    return maxc;
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

void localsearch(int maxStep){
    step=1;
    int i,j,k,h,l,v,c;
    int best_in_c;
    int maxc;
    int flag = 0;
    int rand_n;
    int init_remove = 0;

    for(v = 0; v < vertex_num; v++)
            uncover_vertex_index[v] = -1;//////////////////////////////////

    while(step<=maxStep)
    {
        i = -1;
        /*
            if(step % 100 == 0)
            {
                        times(&finish);
                        double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
                        finish_time = round(finish_time * 100)/100.0;
                        printf(" %ld %ld %f ", current_sol, best_value, finish_time);
            buger();
            }
        */
        while(uncover_num == 0)
        {
            if(t_length > 0)
            {

                    rand_n = rand()%t_length;
                    i = t[rand_n];
            }
            else
            {
                    update_best_sol();///////////////////
                    //i = find_best_in_c_simp(1);
                    //if(state == 0)
                        i = find_best_in_c_simp(0);
                    init_remove = 1;
            }
            remove(i, init_remove);
        }

        times(&finish);
        double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        finish_time = round(finish_time * 100)/100.0;
        if(finish_time>time_limit) break;

        init_remove = 1;
        remove_num1 = i;
        // cout << "r " << i;
        //memset(uncover_vertex_index, -1, sizeof(uncover_vertex_index));
        uncover_length = 0;

        uncover_vertex[uncover_length] = remove_num1;
        uncover_vertex_index[remove_num1] = uncover_length;
        uncover_length++;

        h = vertex_neightbourNum[remove_num1];
        for(v = 0; v < h;v++)
        {
            i = vertex[remove_num1][v];
            if(cs[i].is_in_c == 1)
                continue;
            if(cs[i].is_in_c == 0 && uncover_vertex_index[i] == -1)
            {
                uncover_vertex[uncover_length] = i;
                uncover_vertex_index[i] = uncover_length;
                uncover_length++;
            }
            k = vertex_neightbourNum[i];
            for(j = 0; j < k; j++)
            {
                c = vertex[i][j];
                if(c == remove_num1)
                    continue;
                if(cs[c].is_in_c == 0 && uncover_vertex_index[c] == -1)
                {
                    uncover_vertex[uncover_length] = c;
                    uncover_vertex_index[c] = uncover_length;
                    uncover_length++;
                }
            }
        }

        best_in_c=find_best_in_c(1);
        if(state==0)best_in_c=find_best_in_c(0);
        tabu_list.clear();
        remove(best_in_c, init_remove);
        cs[best_in_c].time_stamp=step;
        //cout << " "<<best_in_c;
        //buger();
        if(uncover_vertex_index[best_in_c] == -1)
        {
            uncover_vertex[uncover_length] = best_in_c;
            uncover_vertex_index[best_in_c] = uncover_length;
            uncover_length++;
        }

        h = vertex_neightbourNum[best_in_c];
        for(v = 0; v < h; v++)
        {
            i = vertex[best_in_c][v];
            if(cs[i].is_in_c == 1)
                    continue;
            if(cs[i].is_in_c == 0 && uncover_vertex_index[i] == -1)
            {
                uncover_vertex[uncover_length] = i;
                uncover_vertex_index[i] = uncover_length;
                uncover_length++;
            }
            k = vertex_neightbourNum[i];
            for(j = 0; j < k; j++)
            {
                c = vertex[i][j];
                if(c == best_in_c)
                        continue;
                if(cs[c].is_in_c == 0 && uncover_vertex_index[c] == -1)
                {
                    uncover_vertex[uncover_length] = c;
                    uncover_vertex_index[c] = uncover_length;
                    uncover_length++;
                }
            }
        }
        //cout << "a";
        while(uncover_num>0)
        {
            int sr=INT_MIN, ct;
            maxc=-1;

            for(v = 0; v < uncover_length; v++)
            {
                j = uncover_vertex[v];
                if(cs[j].config == 0 || cs[j].is_in_c == 1)
                        continue;
                k=compare(sr,ct, cs[j].score, cs[j].cost);
                if(sr==INT_MIN||k<0)
                {
                    sr=cs[j].score;
                    ct=cs[j].cost;
                    maxc=j;
                }
                else if(k == 0 &&  cs[maxc].config == cs[j].config && cs[j].time_stamp < cs[maxc].time_stamp || k == 0 && cs[maxc].config < cs[j].config)
                {
                    maxc=j;
                }
            }
            assert(maxc != -1);
            add(maxc, 1, 1);
            tabu_list.insert(maxc);
            cs[maxc].time_stamp = step;
            uncov_r_weight_inc();
            //cout <<" " <<maxc ;
        }
        //buger();
        for(v = 0; v < uncover_length; v++)
                uncover_vertex_index[uncover_vertex[v]] = -1;

        //update_best_sol();
        step++;
        each_iter++;
        bms_thre = bms_thre*2;
        if(bms_thre >cs_size)
            bms_thre-=cs_size;
        if(step>=total_step) 
            break;
    }
}


void update_best_sol(){
    int i,j,v;
    if(current_sol < best_value)
    {
	    each_iter=0;
        best_value = current_sol;
        for(v=0;v<remain_num;v++)
        {
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
	if(argc<3){
        printf("input wrong\n");
        return 0;
	}

	build_instance_massive(argv[1]);
  	//printf("%s	",argv[1]);
    seed=atoi(argv[2]);
	time_limit=atof(argv[3]);
	total_step=INT_MAX;

	srand(seed);

    init();
	init_reduce();
	/*printf("c %s",argv[1]);
	printf("p cnf %d %d 201\n", remain_num, remain_num*2);
	for(int i = 0; i < remain_num; i++){
    int v1 = remain_vertex[i];
    printf("201 %d ", v1+1);
    for (int j = 0; j < vertex_neightbourNum[v1]; j++){
			int v2 = vertex[v1][j];
			printf("%d ",v2+1);
    }
    printf("0\n");
	}

	for(int i = 0; i < remain_num; i++){
    int v1 = remain_vertex[i];
    printf("1 -%d 0\n",v1+1);

	}
	exit(0);
	*/


	times(&start);
	start_time = start.tms_utime + start.tms_stime;
	bms = 100;

	init_best();

	localsearch(total_step);

	if(!check()) printf("wrong answer \n");
	printf("%ld	",best_value);
	printf("%.2f\n",  real_time);
	free_all();
	return 0;
}
