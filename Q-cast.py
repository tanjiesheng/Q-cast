#    SimQN: a discrete-event simulator for the quantum networks
#    Copyright (C) 2021-2022 Lutong Chen, Jian Li, Kaiping Xue
#    University of Science and Technology of China, USTC.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from typing import Callable, Dict, List, Tuple, Union,Set
import math
from qns.network.requests import Request
from qns.entity.node.node import QNode
from qns.entity.qchannel.qchannel import QuantumChannel
from qns.entity.cchannel.cchannel import ClassicChannel
from qns.network.route.route import RouteImpl, NetworkRouteError
from qns.models.epr import WernerStateEntanglement
from qns.utils.rnd import get_randint, get_rand
import copy
'最近修改了将纠缠资源按保真度降序进行排序'
"添加了按需产生纠缠资源的路由算法"



class MyRouteAlgorithm(RouteImpl):
    """
    This is my route algorithm implement
    """
    INF = math.inf
    def __init__(self, name: str = "My",
                 metric_func: Callable[[Union[QuantumChannel, ClassicChannel]], float] = None) -> None:
        """
        Args:
            name: the routing algorithm's name
            metric_func: the function that returns the metric for each channel.
                The default is the const function m(l)=1
        """
        self.name = name
        self.route_table = {}
        self.quantum_table = {}
        if metric_func is None:
            self.metric_func = lambda _: 1
        else:
            self.metric_func = metric_func


    def build(self, nodes: List[QNode], channels: List[Union[QuantumChannel, ClassicChannel]]):

        for n in nodes:
            selected = []
            unselected = [u for u in nodes]

            d = {}
            for nn in nodes:
                if nn == n:
                    d[n] = [0, []]
                else:
                    d[nn] = [self.INF, [nn]]

            while len(unselected) != 0:
                ms = unselected[0]
                mi = d[ms][0]

                for s in unselected:
                    if d[s][0] < mi:
                        ms = s
                        mi = d[s][0]

                # d[ms] = [d[ms][0], d[ms][1]]
                selected.append(ms)
                unselected.remove(ms)

                for link in channels:
                    if ms not in link.node_list:
                        continue
                    if len(link.node_list) < 2:
                        raise NetworkRouteError("broken link")
                    idx = link.node_list.index(ms)
                    idx_s = 1 - idx
                    s = link.node_list[idx_s]
                    if s in unselected and d[s][0] > d[ms][0] + self.metric_func(link):
                        d[s] = [d[ms][0] + self.metric_func(link), [ms] + d[ms][1]]

            for nn in nodes:
                d[nn][1] = [nn] + d[nn][1]
            self.route_table[n] = d
            #quantum source

    def minlen(self,current:QNode,dest:QNode):
        """
        estimate the fidelity if current node choose the shortest path to build the E2E entanglement

        """
        n1 = self.route_table[current]
        for key in n1:
            if key == dest:
                minlen = n1[key][0]
        return minlen

    def estimate(self,src:QNode,current:QNode,dest:QNode,fidelity: float,epr:WernerStateEntanglement =None):
        """

        current: current node
        dest: destnation node
        return: the most optimal fidelity while using the shortest path to swap entanglement between curr and dest
        """
        minlen = self.minlen(current,dest)
        F = 0.0
        n = 0.0
        m = 0.0
        w= 0.0
        if current==src:
            w1 = (fidelity * 4 - 1) / 3
            w = pow(w1,minlen)
            F = (w * 3 + 1) / 4
            return F

        else:                                                  ##直接用之前保留的epr资源

            w = epr.w
            l = self.minlen(src,current)
            w = w*pow(w,minlen/l)
            F = (w * 3 + 1) / 4
            return F


    def query(self, src: QNode, dest: QNode,lock_table: Dict[Request,WernerStateEntanglement],Fth: float) -> List[Tuple[float, float, List[QNode]]]:
        """
        query the metric, nexthop and the path

        Args:
            src: the source node
            dest: the destination node

        Returns:
            A list of route paths. The result should be sortted by the priority.
            The element is a tuple containing: metric, the fidelity and the whole path.
        """

        # find the candidate of next hop
        current=src
        # print("qroute_table is",current.qroute_table)
        # for k,v in current.qroute_table.items():
        #     if current.qroute_table[k] == []:
        #         return
        path = [current]
        length = 0
        epr: WernerStateEntanglement = None  # 标记延长的纠缠资源
        epr1:WernerStateEntanglement=None
        m = 1
        while(current!=dest or current!=None):
            #qroute_table = copy.deepcopy(current.qroute_table)
            m0 = m
            fidelity = 0.97  ##初始纠缠资源的保真度和初始需要预留的资源数量
            F = self.estimate(src, current, dest, fidelity,epr)  ##预估当前节点到目的节点的保真度
            Fth2 = (Fth ** 2 + (1 - Fth) ** 2 / 9) / \
                       (Fth ** 2 + 5 / 9 * (1 - Fth) ** 2 + 2 / 3 * Fth * (1 - Fth))
            #Fth = Fth2  ##  设置冗余（从Fth到Fth2）
            while (F < Fth):
                m0 = m0 + m0
                F = (F * F + (1 - F) * (1 - F) / 9) / \
                (F * F + 5 / 9 * (1 - F) * (1 - F) + 1 / 3 * F * (1 - F) + 1 / 3 * F * (1 - F))
            if m0 > src.memories[0].capacity:
                m0 = src.memories[0].capacity
            #print("m =",m) ##得到m，m是在寻找下一跳时需要的资源
            if m0 > m:
                m = m0
            nei1 = set()
            candidates:Set(QNode) = set()
            candidate:Set(QNode) =set()
            nei2 = set()
            n1 = self.route_table[current]
            n2 = self.route_table[dest]
            for key in n2:
                nei2.add(key)
            for key in n1:
                if (n1[key][0] == 1):
                    nei1.add(key)
                # find the neighbors of current
            for key in n1:
                #print("error")
                if (key == dest):
                    minlen = n1[key][0]   #源节点和目的节点之间的距离
            C = nei1 & nei2                #current的邻居中能到dest的节点集合
            for key in n2:
                if (key in C and n2[key][0] < minlen):
                    candidate.add(key)
            #print("candidates is ",candidates) ##找到current的候选集,然后从中找到有m个资源的下一跳
            for node in candidate:
                if len(current.qroute_table[node])>=m:
                    candidates.add(node)
            if len(candidates)<=0:
                return #没有满足资源的下一跳
            quantum_source:Dict[QNode,List[WernerStateEntanglement]] = {}
            for i in candidates:
                if i==dest:#目的节点为当前节点的邻居
                    if len(current.qroute_table[dest])>=m:
                        l = len(current.qroute_table[dest])
                        ql:List[WernerStateEntanglement] = current.qroute_table.get(dest,None)
                        quantum_source[dest]  = ql
                         #找到有m份资源的候选节点
                        for k in range(0, len(current.qroute_table[dest])):
                            for j in range(k + 1, len(quantum_source[dest])):
                                if quantum_source[dest][k].w < quantum_source[dest][j].w:
                                    quantum_source[dest][k], quantum_source[dest][j] = quantum_source[dest][j], \
                                                                                       quantum_source[dest][k]
                        #将资源进行保真度排序
                        # next_hop=dest
                        if(m==1):
                            if current == src:
                                epr = quantum_source[dest][0]  #简单的赋值
                                current.memories[0].read(quantum_source[dest][0])
                                dest.memories[0].read(quantum_source[dest][0])
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]  # 选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                            else:
                                # z = quantum_source[next_hop][0]

                                current.memories[0].read(quantum_source[dest][0])
                                dest.memories[0].read(quantum_source[dest][0])
                                epr = epr.swapping(quantum_source[dest][0])

                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    epr = None
                                    #print("fail to swap")
                                    return  # 纠缠交换失败，资源作废，请求失败
                        elif m>1:
                            m1 = m-1

                            for i in range(0,m1):
                                current.memories[0].read(quantum_source[dest][m - i - 1])
                                dest.memories[0].read(quantum_source[dest][m - i - 1])
                                current.memories[0].read(quantum_source[dest][m-i-2])
                                dest.memories[0].read(quantum_source[dest][m-i-2])
                                a = quantum_source[dest][m-i-2]
                                b = quantum_source[dest][m - i - 1]
                                quantum_source[dest][m-i-2] = a.distillation(b)

                                del dest.qroute_table[current][m-i-1]
                                del current.qroute_table[dest][m-i-1]
                                epr1=quantum_source[dest][m-i-2]
                                if quantum_source[dest][m - i - 2] == -1:
                                    del dest.qroute_table[current][m - i - 2]
                                    del current.qroute_table[dest][m - i - 2]
                                        # print("fail to distill")
                                        # print("最后一跳纯化失败1")
                                    break
                                # m = int(m / 2)
                            if epr1 == -1:
                                break  ##此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤

                                # current.memories[0].write(quantum_source[next_hop][0])
                                # next_hop.memories[0].write(quantum_source[next_hop][0])
                            if current == src:
                                epr = epr1
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]
                                # quantum_source[dest].remove(quantum_source[dest][0])
                                for key in lock_table:
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                            else:
                                # print(epr)
                                epr = epr.swapping(epr1)  # epr 指向延长的纠缠资源
                                # quantum_source[dest].remove(quantum_source[dest][0])
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    #print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                                    return
                        current =dest
                        path.append(dest)
                        length = length+1
                        # print("epr is", epr)
                        # print("fidelity is", epr.fidelity)
                        if epr.fidelity < Fth:
                            return
                        nexthop = path[1]
                        print(epr.fidelity)
                        return [(length, nexthop, path)]                  ##此时成功

                    else:
                        return             #当资源不足的时候，request无法满足，后续可以考虑修改为尽力而为

            if epr1==-1:
                # print("最后一跳纯化失败")
                epr1 =None
                continue

            for i in candidates:
                l = len(current.qroute_table[i])
                if l>=m:
                    ql:List[WernerStateEntanglement] = current.qroute_table.get(i,None)
                    quantum_source[i] = ql  # 找到有m份资源的候选节点

                               #候选节点中没有资源，应该结束event，再继续生成纠缠
            # for key1, value1 in quantum_source.items():
            #     for key2, value2 in lock_table.items():
            #         for i in value1:
            #             if i == value2 and i != None:
            #                 value1.remove(i)                  ##将其他请求占用的第一跳纠缠资源去掉
            # for key,value in quantum_source.items():
            #     if len(quantum_source[key])<m:
            #         print("资源数量为",len(quantum_source[key]))
            #         quantum_source.pop(quantum_source[key])                 ##去掉后如果资源数小于m，去除该节点,这里有可能为空
            if quantum_source =={}:
                #print("lack of qsource")
                return
            for item,item2 in quantum_source.items():
                for i in range(0,len(quantum_source[item])):
                    for j in range(i + 1, len(quantum_source[item])):
                        if quantum_source[item][i].w < quantum_source[item][j].w:
                            quantum_source[item][i], quantum_source[item][j] = quantum_source[item][j], \
                                                                                   quantum_source[item][i]
            #print(quantum_source)  ##得到最终的资源结果,对每个节点的纠缠资源按照保真度降序排列
            node_set: Set(QNode) = set()
            if m == 1:
                l = math.inf
                for key,value in quantum_source.items():
                    if l > self.minlen(key, dest):
                        l = self.minlen(key, dest)
                        next_hop = key  ##实际上应该都是一样的，都是最短路径
                for key,value in quantum_source.items():
                    if self.minlen(key, dest) == l:
                        node_set.add(key)
                w1 = 0
                for key,value in quantum_source.items():
                    if not key in node_set:
                        #print("节点选取错误")
                        quantum_source.pop(key)
                for item in node_set:
                    if w1 < quantum_source[item][0].w :
                        w1 = quantum_source[item][0].w
                        next_hop = item  ##从拥有资源的最短距离的节点中找到保真度最高的节点，即为下一跳

                if current == src:
                    epr = quantum_source[next_hop][0]
                    next_hop.qroute_table[current].remove(epr)
                    current.qroute_table[next_hop].remove(epr)  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表

                    current.memories[0].read(epr)
                    next_hop.memories[0].read(epr)
                    for key,value in lock_table.items():
                        if (key.src==src and key.dest == dest):
                            lock_table[key] = epr
                else:
                    z = quantum_source[next_hop][0]

                    current.memories[0].read(z)
                    next_hop.memories[0].read(z)
                    epr = epr.swapping(z)
                    del current.qroute_table[next_hop][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                    del next_hop.qroute_table[current][0]
                    for key,value in lock_table.items():
                        if (key.src == src and key.dest == dest):
                            lock_table[key] = []   # 将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                    if epr == -1:
                        epr = None
                       #print("fail to swap")
                        return  ##纠缠交换失败，资源作废，请求失败
                #print("epr is",epr)
            elif m > 1:
                l = math.inf
                for key,value in quantum_source.items():
                    if l > self.minlen(key, dest):
                        l = self.minlen(key, dest)
                        next_hop = key
                for key,value in quantum_source.items():
                    if self.minlen(key, dest) == l:
                        node_set.add(key)
                for key,value in quantum_source.items():
                    if not key in node_set:
                       #print("节点选取错误")
                        quantum_source.pop(key)
                        return
                F1 = 0
                for item in node_set:
                    fidelity_dic:Dict[QNode,float]={}
                    m1 = m
                    fidelity_dic[item] = []
                    for i in range(m1):
                        fidelity_dic[item].append(quantum_source[item][i].fidelity)

                    m2 = m-1
                    for i in range(0,m2):
                        a = fidelity_dic[item][m-i-1]
                        b = fidelity_dic[item][m-i-2]
                        fidelity_dic[item][m-i-2] = (a * b + (1 - a) * (1 - b) / 9) / \
                                (a * b + 5 / 9 * (1 - a) * (1 - b) + 1 / 3 * a * (1 - b) + 1 / 3 * b * (1 - a))


                    if fidelity_dic[item][0] > F1:  ##每一个item都有一个最终的理论fidelity
                        F1 = fidelity_dic[item][0]
                        next_hop = item  ##从拥有资源的最短距离的节点中找到保真度最高的节点，即为下一跳

                # num = int(m/2)
                m1 = m-1
                for i in range(0,m1):  #range不包括结束的节点
                    current.memories[0].read(quantum_source[next_hop][m-i-1])
                    next_hop.memories[0].read(quantum_source[next_hop][m-i-1])
                    # s1 = quantum_source[next_hop][i]

                    #s = quantum_source[next_hop][num+i]
                    current.memories[0].read(quantum_source[next_hop][m-i-2])
                    next_hop.memories[0].read(quantum_source[next_hop][m-i-2])##将最终选定的节点中需要用到的资源从节点取出，并在qroutetable上删除
                    # current.qroute_table[next_hop].remove(quantum_source[next_hop][m-i-1])
                    # next_hop.qroute_table[current].remove(quantum_source[next_hop][m-i-1])
                    #current.qroute_table[next_hop].remove(quantum_source[next_hop][num+i])
                    #next_hop.qroute_table[current].remove(quantum_source[next_hop][num+i])
                    a = quantum_source[next_hop][m-i-2]
                    b = quantum_source[next_hop][m - i - 1]

                    quantum_source[next_hop][m - i - 2] = a.distillation(b)
                    del next_hop.qroute_table[current][m-i-1]
                    del current.qroute_table[next_hop][m-i-1]
                    epr1 = quantum_source[next_hop][m-i-2]
                    if epr1 == -1:
                        del next_hop.qroute_table[current][m-i-2]
                        del current.qroute_table[next_hop][m-i-2]
                        #print("fail to distill")

                        break
                  # 此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤
                if epr1 ==-1:
                    epr1 = None
                    continue
                    # current.memories[0].write(quantum_source[next_hop][0])
                    # next_hop.memories[0].write(quantum_source[next_hop][0])
                if current == src:
                    epr = epr1
                    del current.qroute_table[next_hop][0]
                    del next_hop.qroute_table[current][0]
                    for key in lock_table:
                        if (key.src==src and key.dest == dest):
                            lock_table[key] = epr
                else:
                    #print(epr)
                    epr = epr.swapping(epr1)  #epr 指向延长的纠缠资源
                    del current.qroute_table[next_hop][0]
                    del next_hop.qroute_table[current][0]
                    for key,value in lock_table.items():
                        if (key.src == src and key.dest == dest):
                            lock_table[key] = []   ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                    if epr == -1:
                        #print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                        return
            if epr.fidelity<Fth:
                return
            current = next_hop
            length = length+1
            path.append(current)
        if path[-1]!=dest:
            return #纠缠交换失败，返回空，等待下一次的请求

        ##最终找到目的节点

class MyRouteAlgorithm1(RouteImpl):
    """
    This is my route algorithm implement
    """
    INF = math.inf
    def __init__(self, name: str = "My",
                 metric_func: Callable[[Union[QuantumChannel, ClassicChannel]], float] = None) -> None:
        """
        Args:
            name: the routing algorithm's name
            metric_func: the function that returns the metric for each channel.
                The default is the const function m(l)=1
        """
        self.name = name
        self.route_table = {}
        self.quantum_table = {}
        if metric_func is None:
            self.metric_func = lambda _: 1
        else:
            self.metric_func = metric_func


    def build(self, nodes: List[QNode], channels: List[Union[QuantumChannel, ClassicChannel]]):

        for n in nodes:
            selected = []
            unselected = [u for u in nodes]

            d = {}
            for nn in nodes:
                if nn == n:
                    d[n] = [0, []]
                else:
                    d[nn] = [self.INF, [nn]]

            while len(unselected) != 0:
                ms = unselected[0]
                mi = d[ms][0]

                for s in unselected:
                    if d[s][0] < mi:
                        ms = s
                        mi = d[s][0]

                # d[ms] = [d[ms][0], d[ms][1]]
                selected.append(ms)
                unselected.remove(ms)

                for link in channels:
                    if ms not in link.node_list:
                        continue
                    if len(link.node_list) < 2:
                        raise NetworkRouteError("broken link")
                    idx = link.node_list.index(ms)
                    idx_s = 1 - idx
                    s = link.node_list[idx_s]
                    if s in unselected and d[s][0] > d[ms][0] + self.metric_func(link):
                        d[s] = [d[ms][0] + self.metric_func(link), [ms] + d[ms][1]]

            for nn in nodes:
                d[nn][1] = [nn] + d[nn][1]
            self.route_table[n] = d
            #quantum source

    def minlen(self,current:QNode,dest:QNode):
        """
        estimate the fidelity if current node choose the shortest path to build the E2E entanglement

        """
        n1 = self.route_table[current]
        for key in n1:
            if key == dest:
                minlen = n1[key][0]
        return minlen

    def estimate(self,src:QNode,current:QNode,dest:QNode,fidelity: float,epr:WernerStateEntanglement =None):
        """

        current: current node
        dest: destnation node
        return: the most optimal fidelity while using the shortest path to swap entanglement between curr and dest
        """
        minlen = self.minlen(current,dest)
        F = 0.0
        n = 0.0
        m = 0.0
        w= 0.0
        if current==src:
            w1 = (fidelity * 4 - 1) / 3
            w = pow(w1,minlen)
            F = (w * 3 + 1) / 4
            return F

        else:                                                  ##直接用之前保留的epr资源

            w = epr.w
            l = self.minlen(src,current)
            w = w*pow(w,minlen/l)
            F = (w * 3 + 1) / 4
            return F


    def query(self, src: QNode, dest: QNode,lock_table: Dict[Request,WernerStateEntanglement],Fth: float) -> List[Tuple[float, float, List[QNode]]]:
        """
        query the metric, nexthop and the path

        Args:
            src: the source node
            dest: the destination node

        Returns:
            A list of route paths. The result should be sortted by the priority.
            The element is a tuple containing: metric, the fidelity and the whole path.
        """

        # find the candidate of next hop
        current=src
        # print("qroute_table is",current.qroute_table)
        path = [current]
        length = 0
        epr: WernerStateEntanglement = None  # 标记延长的纠缠资源
        epr1:WernerStateEntanglement = None
        m = 1
        while(current!=dest and current!=None):
            #qroute_table = copy.deepcopy(current.qroute_table)
            m0 = 1
            fidelity = 0.97  ##初始纠缠资源的保真度和初始需要预留的资源数量
            F = self.estimate(src, current, dest, fidelity,epr)  # 预估当前节点到目的节点的保真度
            Fth2 = (Fth ** 2 + (1 - Fth) ** 2 / 9) / \
                       (Fth ** 2 + 5 / 9 * (1 - Fth) ** 2 + 2 / 3 * Fth * (1 - Fth))
            #Fth = Fth2  ##  设置冗余（从Fth到Fth2）
            while (F < Fth):
                m0 = m0 + m0
                F = (F * F + (1 - F) * (1 - F) / 9) / \
                (F * F + 5 / 9 * (1 - F) * (1 - F) + 1 / 3 * F * (1 - F) + 1 / 3 * F * (1 - F))
            if m0 > src.memories[0].capacity:
                m0 = src.memories[0].capacity
            #print("m =",m) ##得到m，m是在寻找下一跳时需要的资源
            if m0 > m:
                m = m0
            nei1 = set()
            candidates:Set(QNode) = set()
            nei2 = set()
            n1 = self.route_table[current]
            n2 = self.route_table[dest]
            for key in n2:
                nei2.add(key)
            for key in n1:
                if (n1[key][0] == 1):
                    nei1.add(key)
                # find the neighbors of current
            for key in n1:
                #print("error")
                if (key == dest):
                    minlen = n1[key][0]   #源节点和目的节点之间的距离
            C = nei1 & nei2                #current的邻居中能到dest的节点集合
            for key in n2:
                if (key in C and n2[key][0] < minlen):
                    candidates.add(key)
            #print("candidates is ",candidates) # 找到current的候选集,然后从中找到有m个资源的下一跳
            quantum_source:Dict[QNode,List[WernerStateEntanglement]] = {}
            for i in candidates:
                if i==dest:#目的节点为当前节点的邻居
                    current.qroute_table[dest] = []
                    dest.qroute_table[current] = []
                    for k in range(0,m*2):
                        name_id=1
                        k1 = get_rand(0, 1)
                        if k1 < 0.5:  # 随机生成纠缠链接，概率0.5
                            name_id = k + 1
                            e1 = WernerStateEntanglement(name=f"e{name_id}", fidelity=fidelity)  # 生成m份纠缠资源
                            current.memories[0].write(e1)  # 写入
                            dest.memories[0].write(e1)
                            current.qroute_table[dest].append(e1)  # 修改qroutetable，将资源放进去
                            dest.qroute_table[current].append(e1)
                    if len(current.qroute_table[dest])>=m:
                        l = len(current.qroute_table[dest])
                        ql:List[WernerStateEntanglement] = current.qroute_table.get(dest,None)
                        quantum_source[dest]  = ql
                         #找到有m份资源的候选节点
                        for k in range(0, len(current.qroute_table[dest])):
                            for j in range(k + 1, len(quantum_source[dest])):
                                if quantum_source[dest][k].w < quantum_source[dest][j].w:
                                    quantum_source[dest][k], quantum_source[dest][j] = quantum_source[dest][j], \
                                                                                       quantum_source[dest][k]
                        #将资源进行保真度排序
                        # next_hop=dest
                        if(m==1):
                            if current == src:
                                epr = quantum_source[dest][0]  #简单的赋值
                                current.memories[0].read(quantum_source[dest][0])
                                dest.memories[0].read(quantum_source[dest][0])
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]  # 选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                            else:
                                # z = quantum_source[next_hop][0]

                                current.memories[0].read(quantum_source[dest][0])
                                dest.memories[0].read(quantum_source[dest][0])
                                epr = epr.swapping(quantum_source[dest][0])

                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    epr = None
                                    #print("fail to swap")
                                    return  # 纠缠交换失败，资源作废，请求失败
                        elif m>1:
                            m1 = m-1

                            for i in range(0,m1):
                                current.memories[0].read(quantum_source[dest][m - i - 1])
                                dest.memories[0].read(quantum_source[dest][m - i - 1])
                                current.memories[0].read(quantum_source[dest][m-i-2])
                                dest.memories[0].read(quantum_source[dest][m-i-2])
                                a = quantum_source[dest][m-i-2]
                                b = quantum_source[dest][m - i - 1]
                                quantum_source[dest][m-i-2] = a.distillation(b)

                                del dest.qroute_table[current][m-i-1]
                                del current.qroute_table[dest][m-i-1]
                                epr1=quantum_source[dest][m-i-2]
                                if quantum_source[dest][m - i - 2] == -1:
                                    del dest.qroute_table[current][m - i - 2]
                                    del current.qroute_table[dest][m - i - 2]
                                        # print("fail to distill")
                                        # print("最后一跳纯化失败1")
                                    break
                                # m = int(m / 2)
                            if epr1 == -1:
                                break  ##此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤

                                # current.memories[0].write(quantum_source[next_hop][0])
                                # next_hop.memories[0].write(quantum_source[next_hop][0])
                            if current == src:
                                epr = epr1
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]
                                # quantum_source[dest].remove(quantum_source[dest][0])
                                for key in lock_table:
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                            else:
                                # print(epr)
                                epr = epr.swapping(epr1)  # epr 指向延长的纠缠资源
                                # quantum_source[dest].remove(quantum_source[dest][0])
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    #print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                                    return
                        current =dest
                        path.append(dest)
                        length = length+1
                        # print("epr is", epr)
                        # print("fidelity is", epr.fidelity)
                        if epr.fidelity < Fth:
                            return
                        nexthop = path[1]
                        return [(length, nexthop, path)]                  ##此时成功
                    else:
                        return             #当资源不足的时候，request无法满足，后续可以考虑修改为尽力而为

            if epr1==-1:
                # print("最后一跳纯化失败")
                epr1 =None
                continue

            for node in candidates:  # 生成资源
                current.qroute_table[node]=[]
                node.qroute_table[current]=[]
                for k in range(0, m * 2):
                    k1 = get_rand(0, 1)
                    if k1 < 0.5:  # 随机生成纠缠链接，概率0.5
                        name_id = k + 1
                        e1 = WernerStateEntanglement(name=f"e{name_id}", fidelity=fidelity)  # 生成m份纠缠资源
                        current.memories[0].write(e1)  # 写入
                        node.memories[0].write(e1)

                        current.qroute_table[node].append(e1)  # 修改qroutetable，将资源放进去
                        node.qroute_table[current].append(e1)
                l = len(current.qroute_table[node])
                if l>=m:
                    ql:List[WernerStateEntanglement] = current.qroute_table.get(node,None)
                    quantum_source[node] = ql  # 找到有m份资源的候选节点

                               #候选节点中没有资源，应该结束event，再继续生成纠缠
            # for key1, value1 in quantum_source.items():
            #     for key2, value2 in lock_table.items():
            #         for i in value1:
            #             if i == value2 and i != None:
            #                 value1.remove(i)                  ##将其他请求占用的第一跳纠缠资源去掉
            # for key,value in quantum_source.items():
            #     if len(quantum_source[key])<m:
            #         print("资源数量为",len(quantum_source[key]))
            #         quantum_source.pop(quantum_source[key])                 ##去掉后如果资源数小于m，去除该节点,这里有可能为空
            if quantum_source =={}:
                #print("lack of qsource")
                return
            for item,item2 in quantum_source.items():
                for i in range(0,len(quantum_source[item])):
                    for j in range(i + 1, len(quantum_source[item])):
                        if quantum_source[item][i].w < quantum_source[item][j].w:
                            quantum_source[item][i], quantum_source[item][j] = quantum_source[item][j], \
                                                                                   quantum_source[item][i]
            #print(quantum_source)  ##得到最终的资源结果,对每个节点的纠缠资源按照保真度降序排列
            node_set: Set(QNode) = set()
            if m == 1:
                l = math.inf
                for key,value in quantum_source.items():
                    if l > self.minlen(key, dest):
                        l = self.minlen(key, dest)
                        next_hop = key  ##实际上应该都是一样的，都是最短路径
                for key,value in quantum_source.items():
                    if self.minlen(key, dest) == l:
                        node_set.add(key)
                w1 = 0
                for key,value in quantum_source.items():
                    if not key in node_set:
                        #print("节点选取错误")
                        quantum_source.pop(key)
                for item in node_set:
                    if w1 < quantum_source[item][0].w :
                        w1 = quantum_source[item][0].w
                        next_hop = item  ##从拥有资源的最短距离的节点中找到保真度最高的节点，即为下一跳

                if current == src:
                    epr = quantum_source[next_hop][0]
                    next_hop.qroute_table[current].remove(epr)
                    current.qroute_table[next_hop].remove(epr)  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表

                    current.memories[0].read(epr)
                    next_hop.memories[0].read(epr)
                    for key,value in lock_table.items():
                        if (key.src==src and key.dest == dest):
                            lock_table[key] = epr
                else:
                    z = quantum_source[next_hop][0]

                    current.memories[0].read(z)
                    next_hop.memories[0].read(z)
                    epr = epr.swapping(z)
                    del current.qroute_table[next_hop][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                    del next_hop.qroute_table[current][0]
                    for key,value in lock_table.items():
                        if (key.src == src and key.dest == dest):
                            lock_table[key] = []   # 将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                    if epr == -1:
                        epr = None
                       #print("fail to swap")
                        return  ##纠缠交换失败，资源作废，请求失败
                #print("epr is",epr)
            elif m > 1:
                l = math.inf
                for key,value in quantum_source.items():
                    if l > self.minlen(key, dest):
                        l = self.minlen(key, dest)
                        next_hop = key
                for key,value in quantum_source.items():
                    if self.minlen(key, dest) == l:
                        node_set.add(key)
                for key,value in quantum_source.items():
                    if not key in node_set:
                       #print("节点选取错误")
                        quantum_source.pop(key)
                        return
                F1 = 0
                for item in node_set:
                    fidelity_dic:Dict[QNode,float]={}
                    m1 = m
                    fidelity_dic[item] = []
                    for i in range(m1):
                        fidelity_dic[item].append(quantum_source[item][i].fidelity)

                    m2 = m-1
                    for i in range(0,m2):
                        a = fidelity_dic[item][m-i-1]
                        b = fidelity_dic[item][m-i-2]
                        fidelity_dic[item][m-i-2] = (a * b + (1 - a) * (1 - b) / 9) / \
                                (a * b + 5 / 9 * (1 - a) * (1 - b) + 1 / 3 * a * (1 - b) + 1 / 3 * b * (1 - a))


                    if fidelity_dic[item][0] > F1:  ##每一个item都有一个最终的理论fidelity
                        F1 = fidelity_dic[item][0]
                        next_hop = item  ##从拥有资源的最短距离的节点中找到保真度最高的节点，即为下一跳

                # num = int(m/2)
                m1 = m-1
                for i in range(0,m1):  #range不包括结束的节点
                    current.memories[0].read(quantum_source[next_hop][m-i-1])
                    next_hop.memories[0].read(quantum_source[next_hop][m-i-1])
                    # s1 = quantum_source[next_hop][i]

                    #s = quantum_source[next_hop][num+i]
                    current.memories[0].read(quantum_source[next_hop][m-i-2])
                    next_hop.memories[0].read(quantum_source[next_hop][m-i-2])##将最终选定的节点中需要用到的资源从节点取出，并在qroutetable上删除
                    # current.qroute_table[next_hop].remove(quantum_source[next_hop][m-i-1])
                    # next_hop.qroute_table[current].remove(quantum_source[next_hop][m-i-1])
                    #current.qroute_table[next_hop].remove(quantum_source[next_hop][num+i])
                    #next_hop.qroute_table[current].remove(quantum_source[next_hop][num+i])
                    a = quantum_source[next_hop][m-i-2]
                    b = quantum_source[next_hop][m - i - 1]

                    quantum_source[next_hop][m - i - 2] = a.distillation(b)
                    del next_hop.qroute_table[current][m-i-1]
                    del current.qroute_table[next_hop][m-i-1]
                    epr1 = quantum_source[next_hop][m-i-2]
                    if epr1 == -1:
                        del next_hop.qroute_table[current][m-i-2]
                        del current.qroute_table[next_hop][m-i-2]
                        #print("fail to distill")

                        break
                  # 此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤
                if epr1 ==-1:
                    epr1 = None
                    continue
                    # current.memories[0].write(quantum_source[next_hop][0])
                    # next_hop.memories[0].write(quantum_source[next_hop][0])
                if current == src:
                    epr = epr1
                    del current.qroute_table[next_hop][0]
                    del next_hop.qroute_table[current][0]
                    for key in lock_table:
                        if (key.src==src and key.dest == dest):
                            lock_table[key] = epr
                else:
                    #print(epr)
                    epr = epr.swapping(epr1)  #epr 指向延长的纠缠资源
                    del current.qroute_table[next_hop][0]
                    del next_hop.qroute_table[current][0]
                    for key,value in lock_table.items():
                        if (key.src == src and key.dest == dest):
                            lock_table[key] = []   ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                    if epr == -1:
                        #print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                        return
            if epr.fidelity<Fth:
                return
            current = next_hop
            length = length+1
            path.append(current)
        if path[-1]!=dest:
            return #纠缠交换失败，返回空，等待下一次的请求

        ##最终找到目的节点

class PurifyDijkstraRouteAlgorithm(RouteImpl):
    """
    This is my route algorithm implement
    """
    INF = math.inf
    def __init__(self, name: str = "Purify",
                 metric_func: Callable[[Union[QuantumChannel, ClassicChannel]], float] = None) -> None:
        """
        Args:
            name: the routing algorithm's name
            metric_func: the function that returns the metric for each channel.
                The default is the const function m(l)=1
        """
        self.name = name
        self.route_table = {}
        self.quantum_table = {}
        if metric_func is None:
            self.metric_func = lambda _: 1
        else:
            self.metric_func = metric_func


    def build(self, nodes: List[QNode], channels: List[Union[QuantumChannel, ClassicChannel]]):

        for n in nodes:
            selected = []
            unselected = [u for u in nodes]

            d = {}
            for nn in nodes:
                if nn == n:
                    d[n] = [0, []]
                else:
                    d[nn] = [self.INF, [nn]]

            while len(unselected) != 0:
                ms = unselected[0]
                mi = d[ms][0]

                for s in unselected:
                    if d[s][0] < mi:
                        ms = s
                        mi = d[s][0]

                # d[ms] = [d[ms][0], d[ms][1]]
                selected.append(ms)
                unselected.remove(ms)

                for link in channels:
                    if ms not in link.node_list:
                        continue
                    if len(link.node_list) < 2:
                        raise NetworkRouteError("broken link")
                    idx = link.node_list.index(ms)
                    idx_s = 1 - idx
                    s = link.node_list[idx_s]
                    if s in unselected and d[s][0] > d[ms][0] + self.metric_func(link):
                        d[s] = [d[ms][0] + self.metric_func(link), [ms] + d[ms][1]]

            for nn in nodes:
                d[nn][1] = [nn] + d[nn][1]
            self.route_table[n] = d
            #quantum source

    def minlen(self,current:QNode,dest:QNode):
        """
        estimate the fidelity if current node choose the shortest path to build the E2E entanglement

        """
        n1 = self.route_table[current]
        for key in n1:
            if key == dest:
                minlen = n1[key][0]
        return minlen


    def query(self, src: QNode, dest: QNode,lock_table: Dict[Request,WernerStateEntanglement],Fth: float) -> List[Tuple[float, float, List[QNode]]]:
        """
        query the metric, nexthop and the path

        Args:
            src: the source node
            dest: the destination node

        Returns:
            A list of route paths. The result should be sortted by the priority.
            The element is a tuple containing: metric, the fidelity and the whole path.
        """

        # find the candidate of next hop
        current=src
        # print("qroute_table is",current.qroute_table)
        for k,v in current.qroute_table.items():
            if current.qroute_table[k] == []:
                return
        path = [current]
        length = 0
        epr: WernerStateEntanglement = None  # 标记延长的纠缠资源
        epr1: WernerStateEntanglement = None
        m = 1
        while(current!=dest or current!=None):
            #qroute_table = copy.deepcopy(current.qroute_table)
            m0 = 1
            if epr != None and epr == -1:
                m0 = m0 + m0 #当保真度不足时翻倍预留资源
                epr = None
            if m0 > src.memories[0].capacity:
                m0 = src.memories[0].capacity
            #print("m =",m) ##得到m，m是在寻找下一跳时需要的资源
            if m0 > m:
                m = m0
            nei1 = set()
            candidates:Set(QNode) = set()
            nei2 = set()
            n1 = self.route_table[current]
            n2 = self.route_table[dest]
            for key in n2:
                nei2.add(key)
            for key in n1:
                if (n1[key][0] == 1):
                    nei1.add(key)
                # find the neighbors of current
            for key in n1:
                #print("error")
                if (key == dest):
                    minlen = n1[key][0]   #源节点和目的节点之间的距离
            C = nei1 & nei2                #current的邻居中能到dest的节点集合
            for key in n2:
                if (key in C and n2[key][0] < minlen):
                    candidates.add(key)
            #print("candidates is ",candidates) ##找到current的候选集,然后从中找到有m个资源的下一跳
            quantum_source:Dict[QNode,List[WernerStateEntanglement]] = {}
            for i in candidates:
                if i==dest:#目的节点为当前节点的邻居
                    if len(current.qroute_table[dest])>=m:
                        l = len(current.qroute_table[dest])
                        ql:List[WernerStateEntanglement] = current.qroute_table.get(dest,None)
                        quantum_source[dest]  = ql
                         #找到有m份资源的候选节点
                        for k in range(0, len(current.qroute_table[dest])):
                            for j in range(k + 1, len(quantum_source[dest])):
                                if quantum_source[dest][k].w < quantum_source[dest][j].w:
                                    quantum_source[dest][k], quantum_source[dest][j] = quantum_source[dest][j], \
                                                                                       quantum_source[dest][k]
                        # #将资源进行保真度排序
                        # next_hop=dest
                        if(m==1):
                            if current == src:
                                epr = quantum_source[dest][0]  #简单的赋值
                                current.memories[0].read(quantum_source[dest][0])
                                dest.memories[0].read(quantum_source[dest][0])
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]  # 选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                            else:
                                # z = quantum_source[next_hop][0]

                                current.memories[0].read(quantum_source[dest][0])
                                dest.memories[0].read(quantum_source[dest][0])
                                epr = epr.swapping(quantum_source[dest][0])

                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    epr = None
                                    #print("fail to swap")
                                    return  # 纠缠交换失败，资源作废，请求失败
                        elif m>1:
                            m1 = m-1
                            for i in range(0,m1):
                                current.memories[0].read(quantum_source[dest][m - i - 1])
                                dest.memories[0].read(quantum_source[dest][m - i - 1])
                                current.memories[0].read(quantum_source[dest][m-i-2])
                                dest.memories[0].read(quantum_source[dest][m-i-2])
                                a = quantum_source[dest][m-i-2]
                                b = quantum_source[dest][m - i - 1]
                                quantum_source[dest][m-i-2] = a.distillation(b)

                                del dest.qroute_table[current][m-i-1]
                                del current.qroute_table[dest][m-i-1]
                                epr1=quantum_source[dest][m-i-2]

                                if quantum_source[dest][m - i - 2] == -1:
                                    del dest.qroute_table[current][m - i - 2]
                                    del current.qroute_table[dest][m - i - 2]
                                        # print("fail to distill")
                                        # print("最后一跳纯化失败1")
                                    break
                                # m = int(m / 2)
                            if epr1 == -1:
                                break  ##此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤

                                # current.memories[0].write(quantum_source[next_hop][0])
                                # next_hop.memories[0].write(quantum_source[next_hop][0])
                            if current == src:
                                epr = epr1
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]
                                # quantum_source[dest].remove(quantum_source[dest][0])
                                for key in lock_table:
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                            else:
                                # print(epr)
                                epr = epr.swapping(epr1)  # epr 指向延长的纠缠资源
                                # quantum_source[dest].remove(quantum_source[dest][0])
                                del dest.qroute_table[current][0]
                                del current.qroute_table[dest][0]
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    #print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                                    return
                        current =dest
                        path.append(dest)
                        length = length+1
                        # print("epr is", epr)
                        # print("fidelity is", epr.fidelity)
                        if epr.fidelity < Fth:
                            epr = -1
                            break
                        nexthop = path[1]
                        print(epr.fidelity)
                        return [(length, nexthop, path)]                  ##此时成功
                    else:
                        return             #当资源不足的时候，request无法满足，后续可以考虑修改为尽力而为
            if epr == -1 or epr1 == -1:
                if epr1==-1: #最后一跳纯化失败
                    epr1 = None
                    continue
                elif epr == -1:#纯化出来的链接达不到保真度要求，返回重新开始，并加性提高资源预留
                    current = src
                    path = [current]
                    length = 0
                    continue
            for i in candidates:
                l = len(current.qroute_table[i])#3333
                if l>=m:
                    ql:List[WernerStateEntanglement] = current.qroute_table.get(i,None)
                    quantum_source[i] = ql  # 找到有m份资源的候选节点

                               #候选节点中没有资源，应该结束event，再继续生成纠缠
            # for key1, value1 in quantum_source.items():
            #     for key2, value2 in lock_table.items():
            #         for i in value1:
            #             if i == value2 and i != None:
            #                 value1.remove(i)                  ##将其他请求占用的第一跳纠缠资源去掉
            # for key,value in quantum_source.items():
            #     if len(quantum_source[key])<m:
            #         print("资源数量为",len(quantum_source[key]))
            #         quantum_source.pop(quantum_source[key])                 ##去掉后如果资源数小于m，去除该节点,这里有可能为空
            if quantum_source =={}:
                #print("lack of qsource")
                return
            for item,item2 in quantum_source.items():
                for i in range(0,len(quantum_source[item])):
                    for j in range(i + 1, len(quantum_source[item])):
                        if quantum_source[item][i].w < quantum_source[item][j].w:
                            quantum_source[item][i], quantum_source[item][j] = quantum_source[item][j], \
                                                                                   quantum_source[item][i]
            #print(quantum_source)  ##得到最终的资源结果,对每个节点的纠缠资源按照保真度降序排列
            node_set: Set(QNode) = set()
            if m == 1:    #m=1的时候肯定在阈值之上
                l = math.inf
                for key,value in quantum_source.items():
                    if l > self.minlen(key, dest):
                        l = self.minlen(key, dest)
                        next_hop = key  ##实际上应该都是一样的，都是最短路径
                for key,value in quantum_source.items():
                    if self.minlen(key, dest) == l:
                        node_set.add(key)
                w1 = 0
                for key,value in quantum_source.items():
                    if not key in node_set:
                        #print("节点选取错误")
                        quantum_source.pop(key)
                for item in node_set:
                    if w1 < quantum_source[item][0].w :
                        w1 = quantum_source[item][0].w
                        next_hop = item  ##从拥有资源的最短距离的节点中找到保真度最高的节点，即为下一跳

                if current == src:
                    epr = quantum_source[next_hop][0]
                    next_hop.qroute_table[current].remove(epr)
                    current.qroute_table[next_hop].remove(epr)  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表

                    current.memories[0].read(epr)
                    next_hop.memories[0].read(epr)
                    for key,value in lock_table.items():
                        if (key.src==src and key.dest == dest):
                            lock_table[key] = epr
                else:
                    z = quantum_source[next_hop][0]

                    current.memories[0].read(z)
                    next_hop.memories[0].read(z)
                    epr = epr.swapping(z)
                    del current.qroute_table[next_hop][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                    del next_hop.qroute_table[current][0]
                    for key,value in lock_table.items():
                        if (key.src == src and key.dest == dest):
                            lock_table[key] = []   # 将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                    if epr == -1:
                        epr = None
                       #print("fail to swap")
                        return  ##纠缠交换失败，资源作废，请求失败
                #print("epr is",epr)
            elif m > 1:
                l = math.inf
                for key,value in quantum_source.items():
                    if l > self.minlen(key, dest):
                        l = self.minlen(key, dest)
                        next_hop = key
                for key,value in quantum_source.items():
                    if self.minlen(key, dest) == l:
                        node_set.add(key)
                for key,value in quantum_source.items():
                    if not key in node_set:
                       #print("节点选取错误")
                        quantum_source.pop(key)
                        return
                F1 = 0
                for item in node_set:
                    fidelity_dic:Dict[QNode,float]={}
                    m1 = m
                    fidelity_dic[item] = []
                    for i in range(m1):
                        fidelity_dic[item].append(quantum_source[item][i].fidelity)

                    m2 = m-1
                    for i in range(0,m2):
                        a = fidelity_dic[item][m-i-1]
                        b = fidelity_dic[item][m-i-2]
                        fidelity_dic[item][m-i-2] = (a * b + (1 - a) * (1 - b) / 9) / \
                                (a * b + 5 / 9 * (1 - a) * (1 - b) + 1 / 3 * a * (1 - b) + 1 / 3 * b * (1 - a))


                    if fidelity_dic[item][0] > F1:  ##每一个item都有一个最终的理论fidelity
                        F1 = fidelity_dic[item][0]
                        next_hop = item  ##从拥有资源的最短距离的节点中找到保真度最高的节点，即为下一跳

                # num = int(m/2)
                m1 = m-1
                for i in range(0,m1):  #range不包括结束的节点
                    current.memories[0].read(quantum_source[next_hop][m-i-1])
                    next_hop.memories[0].read(quantum_source[next_hop][m-i-1])
                    # s1 = quantum_source[next_hop][i]

                    #s = quantum_source[next_hop][num+i]
                    current.memories[0].read(quantum_source[next_hop][m-i-2])
                    next_hop.memories[0].read(quantum_source[next_hop][m-i-2])##将最终选定的节点中需要用到的资源从节点取出，并在qroutetable上删除
                    # current.qroute_table[next_hop].remove(quantum_source[next_hop][m-i-1])
                    # next_hop.qroute_table[current].remove(quantum_source[next_hop][m-i-1])
                    #current.qroute_table[next_hop].remove(quantum_source[next_hop][num+i])
                    #next_hop.qroute_table[current].remove(quantum_source[next_hop][num+i])
                    a = quantum_source[next_hop][m-i-2]
                    b = quantum_source[next_hop][m - i - 1]

                    quantum_source[next_hop][m - i - 2] = a.distillation(b)
                    del next_hop.qroute_table[current][m-i-1]
                    del current.qroute_table[next_hop][m-i-1]
                    epr1 = quantum_source[next_hop][m-i-2]
                    if epr1 == -1:
                        del next_hop.qroute_table[current][m-i-2]
                        del current.qroute_table[next_hop][m-i-2]
                        #print("fail to distill")

                        break
                  # 此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤
                if epr1 ==-1:
                    epr1 = None
                    continue
                    # current.memories[0].write(quantum_source[next_hop][0])
                    # next_hop.memories[0].write(quantum_source[next_hop][0])
                if current == src:
                    epr = epr1
                    del current.qroute_table[next_hop][0]
                    del next_hop.qroute_table[current][0]
                    for key in lock_table:
                        if (key.src==src and key.dest == dest):
                            lock_table[key] = epr
                else:
                    #print(epr)
                    epr = epr.swapping(epr1)  #epr 指向延长的纠缠资源
                    del current.qroute_table[next_hop][0]
                    del next_hop.qroute_table[current][0]
                    for key,value in lock_table.items():
                        if (key.src == src and key.dest == dest):
                            lock_table[key] = []   ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                    if epr == -1:
                        #print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                        return
            if epr.fidelity<Fth:
                epr = -1
                current = src
                length = 0
                path = [current]
                continue
            current = next_hop
            length = length+1
            path.append(current)
        if path[-1]!=dest:
            return #纠缠交换失败，返回空，等待下一次的请求

        ##最终找到目的节点

class DijkstraRouteAlgorithm(RouteImpl):
    """
    This is the qcast route algorithm implement(extended dijkstra with recovery path)
    """

    INF = math.inf

    def __init__(self, name: str = "dijkstra",
                 metric_func: Callable[[Union[QuantumChannel, ClassicChannel]], float] = None) -> None:
        """
        Args:
            name: the routing algorithm's name
            metric_func: the function that returns the metric for each channel.
                The default is the const function m(l)=1
        """
        self.name = name
        self.route_table = {}
        if metric_func is None:
            self.metric_func = lambda _: 1
        else:
            self.metric_func = metric_func

    def build(self, nodes: List[QNode], channels: List[Union[QuantumChannel, ClassicChannel]]):

        for n in nodes:
            selected = []
            unselected = [u for u in nodes]

            d = {}
            for nn in nodes:
                if nn == n:
                    d[n] = [0, []]
                else:
                    d[nn] = [self.INF, [nn]]

            while len(unselected) != 0:
                ms = unselected[0]
                mi = d[ms][0]

                for s in unselected:
                    if d[s][0] < mi:
                        ms = s
                        mi = d[s][0]

                # d[ms] = [d[ms][0], d[ms][1]]
                selected.append(ms)
                unselected.remove(ms)

                for link in channels:
                    if ms not in link.node_list:
                        continue
                    if len(link.node_list) < 2:
                        raise NetworkRouteError("broken link")
                    idx = link.node_list.index(ms)
                    idx_s = 1 - idx
                    s = link.node_list[idx_s]
                    if s in unselected and d[s][0] > d[ms][0] + self.metric_func(link):
                        d[s] = [d[ms][0] + self.metric_func(link), [ms] + d[ms][1]]

            for nn in nodes:
                d[nn][1] = [nn] + d[nn][1]
            self.route_table[n] = d

    def minlen(self,current:QNode,dest:QNode):
        """
        estimate the fidelity if current node choose the shortest path to build the E2E entanglement

        """
        n1 = self.route_table[current]
        for key in n1:
            if key == dest:
                minlen = n1[key][0]
        return minlen

    def query(self, src: QNode, dest: QNode,lock_table: Dict[Request,WernerStateEntanglement],Fth:float) -> List[Tuple[float, QNode, List[QNode]]]:
        """
        query the metric, nexthop and the path

        Args:
            src: the source node
            dest: the destination node

        Returns:
            A list of route paths. The result should be sortted by the priority.
            The element is a tuple containing: metric, the next-hop and the whole path.
        """
        # find the candidate of next hop
        current = src
        # print("qroute_table is",current.qroute_table)
        for k, v in current.qroute_table.items():
            if current.qroute_table[k] == []:
                return
        path = [current]
        length = 0
        epr: WernerStateEntanglement = None  # 标记延长的纠缠资源
        while (current != dest or current != None):
            m=1
            # print("m =",m) ##得到m，m是在寻找下一跳时需要的资源
            nei1 = set()
            candidates: Set(QNode) = set()
            nei2 = set()
            n1 = self.route_table[current]
            n2 = self.route_table[dest]
            for key in n2:
                nei2.add(key)
            for key in n1:
                if (n1[key][0] == 1):
                    nei1.add(key)
                # find the neighbors of current
            for key in n1:
                # print("error")
                if (key == dest):
                    minlen = n1[key][0]  # 源节点和目的节点之间的距离
            C = nei1 & nei2  # current的邻居中能到dest的节点集合
            for key in n2:
                if (key in C and n2[key][0] < minlen):
                    candidates.add(key)
            # print("candidates is ",candidates) ##找到current的候选集,然后从中找到有m个资源的下一跳
            quantum_source: Dict[QNode, List[WernerStateEntanglement]] = {}
            for i in candidates:
                if i == dest:  # 目的节点为当前节点的邻居
                    if len(current.qroute_table[dest]) >= m:
                        l = len(current.qroute_table[dest])
                        ql: List[WernerStateEntanglement] = current.qroute_table.get(dest, None)
                        quantum_source[dest] = ql
                        # 找到有m份资源的候选节点
                        for k in range(0, len(current.qroute_table[dest])):
                            for j in range(k + 1, len(quantum_source[dest])):
                                if quantum_source[dest][k].w < quantum_source[dest][j].w:
                                    quantum_source[dest][k], quantum_source[dest][j] = quantum_source[dest][j], \
                                                                                       quantum_source[dest][k]
                        # #将资源进行保真度排序
                        # next_hop=dest

                        if current == src:
                            epr = quantum_source[dest][0]  # 简单的赋值
                            current.memories[0].read(quantum_source[dest][0])
                            dest.memories[0].read(quantum_source[dest][0])
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = epr
                            del dest.qroute_table[current][0]
                            del current.qroute_table[dest][0]  # 选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                        else:
                                # z = quantum_source[next_hop][0]

                            current.memories[0].read(quantum_source[dest][0])
                            dest.memories[0].read(quantum_source[dest][0])
                            epr = epr.swapping(quantum_source[dest][0])

                            del dest.qroute_table[current][0]
                            del current.qroute_table[dest][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                            if epr == -1:
                                epr = None
                                # print("fail to swap")
                                return  # 纠缠交换失败，资源作废，请求失败
                        current = dest
                        path.append(dest)
                        length = length + 1
                        # print("epr is", epr)
                        # print("fidelity is", epr.fidelity)
                        if epr.fidelity < Fth:
                            return
                        nexthop = path[1]
                        return [(length, nexthop, path)]  ##此时成功
                    else:
                        return  # 当资源不足的时候，request无法满足，后续可以考虑修改为尽力而为
            for i in candidates:
                l = len(current.qroute_table[i])
                if l >= m:
                    ql: List[WernerStateEntanglement] = current.qroute_table.get(i, None)
                    quantum_source[i] = ql  # 找到有m份资源的候选节点

                    # 候选节点中没有资源，应该结束event，再继续生成纠缠
            # for key1, value1 in quantum_source.items():
            #     for key2, value2 in lock_table.items():
            #         for i in value1:
            #             if i == value2 and i != None:
            #                 value1.remove(i)                  ##将其他请求占用的第一跳纠缠资源去掉
            # for key,value in quantum_source.items():
            #     if len(quantum_source[key])<m:
            #         print("资源数量为",len(quantum_source[key]))
            #         quantum_source.pop(quantum_source[key])                 ##去掉后如果资源数小于m，去除该节点,这里有可能为空
            if quantum_source == {}:
                #print("lack of qsource")
                return
            for item, item2 in quantum_source.items():
                for i in range(0, len(quantum_source[item])):
                    for j in range(i + 1, len(quantum_source[item])):
                        if quantum_source[item][i].w < quantum_source[item][j].w:
                            quantum_source[item][i], quantum_source[item][j] = quantum_source[item][j], \
                                                                               quantum_source[item][i]
            # print(quantum_source)  ##得到最终的资源结果,对每个节点的纠缠资源按照保真度降序排列
            node_set: Set(QNode) = set()
            l = math.inf
            for key, value in quantum_source.items():
                if l > self.minlen(key, dest):
                    l = self.minlen(key, dest)
                    next_hop = key  ##实际上应该都是一样的，都是最短路径
            for key, value in quantum_source.items():
                if self.minlen(key, dest) == l:
                    node_set.add(key)
            w1 = 0
            for key, value in quantum_source.items():
                if not key in node_set:
                    print("节点选取错误")
                    quantum_source.pop(key)
            for item in node_set:
                if w1 < quantum_source[item][0].w:
                    w1 = quantum_source[item][0].w
                    next_hop = item  ##从拥有资源的最短距离的节点中找到保真度最高的节点，即为下一跳

            if current == src:
                epr = quantum_source[next_hop][0]
                current.qroute_table[next_hop].remove(epr)  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                next_hop.qroute_table[current].remove(epr)
                current.memories[0].read(epr)
                next_hop.memories[0].read(epr)
                for key, value in lock_table.items():
                    if (key.src == src and key.dest == dest):
                        lock_table[key] = epr
            else:
                z = quantum_source[next_hop][0]

                current.memories[0].read(z)
                next_hop.memories[0].read(z)
                epr = epr.swapping(z)
                current.qroute_table[next_hop].remove(z)  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                next_hop.qroute_table[current].remove(z)
                for key, value in lock_table.items():
                    if (key.src == src and key.dest == dest):
                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                if epr == -1:
                    epr = None
                    # print("fail to swap")
                    return  ##纠缠交换失败，资源作废，请求失败
            if epr.fidelity<Fth:
                return
            current = next_hop
            length = length + 1
            path.append(current)


class QCASTRouteAlgorithm(RouteImpl):
    """
    This is the dijkstra route algorithm implement
    """

    INF = math.inf

    def __init__(self, name: str = "qcast",
                 metric_func: Callable[[Union[QuantumChannel, ClassicChannel]], float] = None) -> None:
        """
        Args:
            name: the routing algorithm's name
            metric_func: the function that returns the metric for each channel.
                The default is the const function m(l)=1
        """
        self.name = name
        self.route_table = {}
        if metric_func is None:
            self.metric_func = lambda _: 1
        else:
            self.metric_func = metric_func

    def build(self, nodes: List[QNode], channels: List[Union[QuantumChannel, ClassicChannel]]):

        for n in nodes:
            selected = []
            unselected = [u for u in nodes]

            d = {}
            for nn in nodes:
                if nn == n:
                    d[n] = [0, []]
                else:
                    d[nn] = [self.INF, [nn]]

            while len(unselected) != 0:
                ms = unselected[0]
                mi = d[ms][0]

                for s in unselected:
                    if d[s][0] < mi:
                        ms = s
                        mi = d[s][0]

                # d[ms] = [d[ms][0], d[ms][1]]
                selected.append(ms)
                unselected.remove(ms)

                for link in channels:
                    if ms not in link.node_list:
                        continue
                    if len(link.node_list) < 2:
                        raise NetworkRouteError("broken link")
                    idx = link.node_list.index(ms)
                    idx_s = 1 - idx
                    s = link.node_list[idx_s]
                    if s in unselected and d[s][0] > d[ms][0] + self.metric_func(link):
                        d[s] = [d[ms][0] + self.metric_func(link), [ms] + d[ms][1]]

            for nn in nodes:
                d[nn][1] = [nn] + d[nn][1]
            self.route_table[n] = d

    def minlen(self,current:QNode,dest:QNode):
        """
        find the length of shortest path to build the E2E entanglement

        """
        n1 = self.route_table[current]
        for key in n1:
            if key == dest:
                minlen = n1[key][0]
        return minlen

    def query(self, src: QNode, dest: QNode,lock_table: Dict[Request,WernerStateEntanglement],Fth:float) -> List[Tuple[float, QNode, List[QNode]]]:
        """
        query the metric, nexthop and the path

        Args:
            src: the source node
            dest: the destination node

        Returns:
            A list of route paths. The result should be sortted by the priority.
            The element is a tuple containing: metric, the next-hop and the whole path.
        """
        current = src
        # find the candidate of next hop
        RecoveryPath = {} #sd元组对key的恢复路径集合
        # print("qroute_table is",current.qroute_table)
        for k, v in current.qroute_table.items():
            if current.qroute_table[k] == []:
                return
        epr: WernerStateEntanglement = None  # 标记延长的纠缠资源
        m = 1

            # print("m =",m) ##得到m，m是在寻找下一跳时需要的资源
        ls: Dict[QNode, List[float, List[QNode]]] = self.route_table.get(src, None)
        if ls is None:
            return
        le = ls.get(dest, None)
        if le is None:
            return
        path: List[QNode] = le[1]
        path = path.copy()
        path.reverse() #主路径path，接下来找恢复路径

        l = len(path)-1
        for i in range(0,l):
            if i+1<=l:
                for j in range(i+1,l+1):
                    nei1 = set()
                    nei2 = set()
                    current = path[i]
                    dest1 = path[j]
                    n1 = self.route_table.get(current,None)
                    n2 = self.route_table.get(dest1,None)
                    for key in n2:
                        nei2.add(key)
                    for key in n1:
                        if (n1[key][0] == 1):
                            nei1.add(key)
                            # find the neighbors of current
                    for key in n1:
                            # print("error")
                        if (key == dest1):
                            minlen = n1[key][0]  # 源节点和目的节点之间的距离
                    C = nei1 & nei2  # current的邻居中能到dest1的节点集合
                    for z in path:
                        C.discard(z) #把路径上的节点去掉
                    sd: Tuple(QNode, QNode) = (current, dest1)
                    Recovery = []
                    for k in C:
                        rpath = n2.get(k,None)
                        if rpath[0]<3:
                            rpath1:List[QNode] = rpath[1]
                            rpath1 = rpath1.copy()
                            rpath1.insert(0,current)
                            Recovery.append(rpath1)

                    RecoveryPath[sd]=Recovery
                        #找到恢复路径，开始建立纠缠链接,先预留资源，再进行纠缠交换
        while (True):
            epr == None #初始化epr
            #if (epr == -1): #需要更多的纠缠资源
            l = len(path)-1
            quantum_source: Dict[QNode, List[WernerStateEntanglement]] = {} #初始化资源情况
            for i in range(0,l): #开始纠缠分发
                j = i+1
                current = path[i]
                next_hop = path[j]
                l1 = len(current.qroute_table[next_hop]) #找纠缠链路
                if l1>=m:
                    ql: List[WernerStateEntanglement]=current.qroute_table.get(next_hop,None)
                    quantum_source[next_hop] = ql  # 找到有m份资源的候选节点
                    for item, item2 in quantum_source.items():
                        for i in range(0, len(quantum_source[item])):
                            for j in range(i + 1, len(quantum_source[item])):
                                if quantum_source[item][i].w < quantum_source[item][j].w:
                                    quantum_source[item][i], quantum_source[item][j] = quantum_source[item][j], \
                                                                                       quantum_source[item][i]
                    # print(quantum_source)  ##得到最终的资源结果,对每个节点的纠缠资源按照保真度降序排列
                    if m == 1:
                        if current == src:
                            epr = quantum_source[next_hop][0]
                            next_hop.qroute_table[current].remove(epr)
                            current.qroute_table[next_hop].remove(epr)  #选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表

                            current.memories[0].read(epr)
                            next_hop.memories[0].read(epr)
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = epr
                        else:
                            z = quantum_source[next_hop][0]
                            current.memories[0].read(z)
                            next_hop.memories[0].read(z)
                            epr = epr.swapping(z)
                            del current.qroute_table[next_hop][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                            del next_hop.qroute_table[current][0]
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = []  # 将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                            if epr == -1:
                                epr = None
                                    # print("fail to swap")
                                return  ##纠缠交换失败，资源作废，请求失败
                            # print("epr is",epr)
                    elif m > 1:
                        while (len(current.qroute_table[next_hop]) >= m):
                            m1 = m - 1
                            for i in range(0, m1):  # range不包括结束的节点
                                current.memories[0].read(quantum_source[next_hop][m - i - 1])
                                next_hop.memories[0].read(quantum_source[next_hop][m - i - 1])
                                # s1 = quantum_source[next_hop][i]

                                # s = quantum_source[next_hop][num+i]
                                current.memories[0].read(quantum_source[next_hop][m - i - 2])
                                next_hop.memories[0].read(
                                    quantum_source[next_hop][m - i - 2])  ##将最终选定的节点中需要用到的资源从节点取出，并在qroutetable上删除
                                # current.qroute_table[next_hop].remove(quantum_source[next_hop][m-i-1])
                                # next_hop.qroute_table[current].remove(quantum_source[next_hop][m-i-1])
                                # current.qroute_table[next_hop].remove(quantum_source[next_hop][num+i])
                                # next_hop.qroute_table[current].remove(quantum_source[next_hop][num+i])
                                a = quantum_source[next_hop][m - i - 2]
                                b = quantum_source[next_hop][m - i - 1]

                                quantum_source[next_hop][m - i - 2] = a.distillation(b)
                                del next_hop.qroute_table[current][m - i - 1]
                                del current.qroute_table[next_hop][m - i - 1]
                                epr1 = quantum_source[next_hop][m - i - 2]
                                if epr1 == -1:
                                    del next_hop.qroute_table[current][m - i - 2]
                                    del current.qroute_table[next_hop][m - i - 2]
                                    # print("fail to distill")
                                    break
                            # 此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤
                            if epr1 == -1:
                                epr1 = None
                                if (len(current.qroute_table[next_hop]) < m):
                                    return
                                continue #纯化失败，重新找资源进行纯化
                                # current.memories[0].write(quantum_source[next_hop][0])
                                # next_hop.memories[0].write(quantum_source[next_hop][0])
                            if current == src:
                                epr = epr1
                                del current.qroute_table[next_hop][0]
                                del next_hop.qroute_table[current][0]
                                for key in lock_table:
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                            else:
                                # print(epr)
                                epr = epr.swapping(epr1)  # epr 指向延长的纠缠资源
                                del current.qroute_table[next_hop][0]
                                del next_hop.qroute_table[current][0]
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    # print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                                    return
                            if epr.fidelity>Fth:
                                break
                            elif (epr.fidelity < Fth):
                                break

                elif l1<m: #资源不足，找恢复路径
                    sd=(current,next_hop)
                    l3 = math.inf
                    rpath4 = []
                    pathlist= RecoveryPath.get(sd, None)
                    if pathlist==None:
                        return #无恢复路径
                    for rpath3 in pathlist:
                        if len(rpath3)<l3 and len(current.qroute_table[rpath3[1]])>=m:
                            l3 = len(rpath3)
                            rpath4 = rpath3
                            #找有足够资源且最短的一条作为恢复路径
                    if rpath4 ==[]:
                        return        #找不到恢复路径，且主路径没资源，请求失败
                    for node in rpath4:
                        path.insert(j,node)#将恢复路径插入原路径
                    l= len(path)-1  #此时路径长度变化了
                    current = path[i]
                    next_hop = path[j]
                    l2 = len(current.qroute_table[next_hop])
                    if l2>=m:
                        ql:List[WernerStateEntanglement]=current.qroute_table.get(next_hop,None)
                        quantum_source[next_hop] = ql  # 找到有m份资源的候选节点
                        for item, item2 in quantum_source.items():
                            for i in range(0, len(quantum_source[item])):
                                for j in range(i + 1, len(quantum_source[item])):
                                    if quantum_source[item][i].w < quantum_source[item][j].w:
                                        quantum_source[item][i], quantum_source[item][j] = quantum_source[item][j], \
                                                                                           quantum_source[item][i]
                        # print(quantum_source)  ##得到最终的资源结果,对每个节点的纠缠资源按照保真度降序排列
                        if m == 1:
                            if current == src:
                                epr = quantum_source[next_hop][0]
                                next_hop.qroute_table[current].remove(epr)
                                current.qroute_table[next_hop].remove(epr)  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                                current.memories[0].read(epr)
                                next_hop.memories[0].read(epr)
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                            else:
                                z = quantum_source[next_hop][0]
                                current.memories[0].read(z)
                                next_hop.memories[0].read(z)
                                epr = epr.swapping(z)
                                del current.qroute_table[next_hop][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                                del next_hop.qroute_table[current][0]
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                            lock_table[key] = []  # 将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    epr = None
                                        # print("fail to swap")
                                    return  ##纠缠交换失败，资源作废，请求失败
                                # print("epr is",epr)
                        elif m > 1:
                            while (len(current.qroute_table[next_hop]) >= m):
                                m1 = m - 1
                                for i in range(0, m1):  # range不包括结束的节点
                                    current.memories[0].read(quantum_source[next_hop][m - i - 1])
                                    next_hop.memories[0].read(quantum_source[next_hop][m - i - 1])
                                    # s1 = quantum_source[next_hop][i]

                                    # s = quantum_source[next_hop][num+i]
                                    current.memories[0].read(quantum_source[next_hop][m - i - 2])
                                    next_hop.memories[0].read(
                                        quantum_source[next_hop][m - i - 2])  ##将最终选定的节点中需要用到的资源从节点取出，并在qroutetable上删除
                                    # current.qroute_table[next_hop].remove(quantum_source[next_hop][m-i-1])
                                    # next_hop.qroute_table[current].remove(quantum_source[next_hop][m-i-1])
                                    # current.qroute_table[next_hop].remove(quantum_source[next_hop][num+i])
                                    # next_hop.qroute_table[current].remove(quantum_source[next_hop][num+i])
                                    a = quantum_source[next_hop][m - i - 2]
                                    b = quantum_source[next_hop][m - i - 1]

                                    quantum_source[next_hop][m - i - 2] = a.distillation(b)
                                    del next_hop.qroute_table[current][m - i - 1]
                                    del current.qroute_table[next_hop][m - i - 1]
                                    epr1 = quantum_source[next_hop][m - i - 2]
                                    if epr1 == -1:
                                        del next_hop.qroute_table[current][m - i - 2]
                                        del current.qroute_table[next_hop][m - i - 2]
                                        # print("fail to distill")
                                        break


                                # 此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤
                                if epr1 == -1:
                                    epr1 = None
                                    if (len(current.qroute_table[next_hop]) < m):  # 资源不足，无法满足请求
                                        return
                                    continue
                                    # current.memories[0].write(quantum_source[next_hop][0])
                                    # next_hop.memories[0].write(quantum_source[next_hop][0])
                                if current == src:
                                    epr = epr1
                                    del current.qroute_table[next_hop][0]
                                    del next_hop.qroute_table[current][0]
                                    for key in lock_table:
                                        if (key.src == src and key.dest == dest):
                                            lock_table[key] = epr
                                else:
                                    # print(epr)
                                    epr = epr.swapping(epr1)  # epr 指向延长的纠缠资源
                                    del current.qroute_table[next_hop][0]
                                    del next_hop.qroute_table[current][0]
                                    for key, value in lock_table.items():
                                        if (key.src == src and key.dest == dest):
                                            lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                    if epr == -1:
                                        # print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                                        return
                                if epr.fidelity > Fth:
                                    break
                                elif (epr.fidelity < Fth):
                                    break
                if epr.fidelity < Fth:
                    m=m+m
                    break
            if(epr.fidelity>Fth):
                return True

class QCASTRouteAlgorithm1(RouteImpl):
    """
    This is the dijkstra route algorithm implement
    """

    INF = math.inf

    def __init__(self, name: str = "qcast",
                 metric_func: Callable[[Union[QuantumChannel, ClassicChannel]], float] = None) -> None:
        """
        Args:
            name: the routing algorithm's name
            metric_func: the function that returns the metric for each channel.
                The default is the const function m(l)=1
        """
        self.name = name
        self.route_table = {}
        if metric_func is None:
            self.metric_func = lambda _: 1
        else:
            self.metric_func = metric_func

    def build(self, nodes: List[QNode], channels: List[Union[QuantumChannel, ClassicChannel]]):

        for n in nodes:
            selected = []
            unselected = [u for u in nodes]

            d = {}
            for nn in nodes:
                if nn == n:
                    d[n] = [0, []]
                else:
                    d[nn] = [self.INF, [nn]]

            while len(unselected) != 0:
                ms = unselected[0]
                mi = d[ms][0]

                for s in unselected:
                    if d[s][0] < mi:
                        ms = s
                        mi = d[s][0]

                # d[ms] = [d[ms][0], d[ms][1]]
                selected.append(ms)
                unselected.remove(ms)

                for link in channels:
                    if ms not in link.node_list:
                        continue
                    if len(link.node_list) < 2:
                        raise NetworkRouteError("broken link")
                    idx = link.node_list.index(ms)
                    idx_s = 1 - idx
                    s = link.node_list[idx_s]
                    if s in unselected and d[s][0] > d[ms][0] + self.metric_func(link):
                        d[s] = [d[ms][0] + self.metric_func(link), [ms] + d[ms][1]]

            for nn in nodes:
                d[nn][1] = [nn] + d[nn][1]
            self.route_table[n] = d

    def minlen(self,current:QNode,dest:QNode):
        """
        find the length of shortest path to build the E2E entanglement

        """
        n1 = self.route_table[current]
        for key in n1:
            if key == dest:
                minlen = n1[key][0]
        return minlen

    def query(self, src: QNode, dest: QNode,lock_table: Dict[Request,WernerStateEntanglement],Fth:float) -> List[Tuple[float, QNode, List[QNode]]]:
        """
        query the metric, nexthop and the path

        Args:
            src: the source node
            dest: the destination node

        Returns:
            A list of route paths. The result should be sortted by the priority.
            The element is a tuple containing: metric, the next-hop and the whole path.
        """
        current = src
        # find the candidate of next hop
        RecoveryPath = {} #sd元组对key的恢复路径集合
        # print("qroute_table is",current.qroute_table)
        epr: WernerStateEntanglement = None  # 标记延长的纠缠资源
        m = 1

            # print("m =",m) ##得到m，m是在寻找下一跳时需要的资源
        ls: Dict[QNode, List[float, List[QNode]]] = self.route_table.get(src, None)
        if ls is None:
            return
        le = ls.get(dest, None)
        if le is None:
            return
        path: List[QNode] = le[1]
        path = path.copy()
        path.reverse() #主路径path，接下来找恢复路径
        path1=path.copy()
        l = len(path)-1
        for i in range(0,l):
            if i+1<=l:
                for j in range(i+1,l+1):
                    nei1 = set()
                    nei2 = set()
                    current = path[i]
                    dest1 = path[j]
                    n1 = self.route_table.get(current,None)
                    n2 = self.route_table.get(dest1,None)
                    for key in n2:
                        nei2.add(key)
                    for key in n1:
                        if (n1[key][0] == 1):
                            nei1.add(key)
                            # find the neighbors of current
                    for key in n1:
                            # print("error")
                        if (key == dest1):
                            minlen = n1[key][0]  # 源节点和目的节点之间的距离
                    C = nei1 & nei2  # current的邻居中能到dest1的节点集合
                    for z in path:
                        C.discard(z) #把路径上的节点去掉
                    sd: Tuple(QNode, QNode) = (current, dest1)
                    Recovery = []
                    for k in C:
                        pathid = 0
                        rpath = n2.get(k, None)  # 点k到点dest1的最短路径
                        for node in rpath[1]:
                            if node == current:
                                pathid = 1
                                break
                        if pathid == 1:
                            pathid = 0
                            continue
                        if rpath[0]<3:
                            rpath1:List[QNode] = rpath[1]
                            rpath1 = rpath1.copy()
                            rpath1.insert(0,current)
                            Recovery.append(rpath1)

                    RecoveryPath[sd]=Recovery
                        #找到恢复路径，开始建立纠缠链接,先预留资源，再进行纠缠交换
        while (m<16):
            i1=0
            epr = None #初始化epr
            #if (epr == -1): #需要更多的纠缠资源
            path=path1.copy()
            current=path[i1]
            next_hop=path[i1]
            l = len(path)-1
            quantum_source: Dict[QNode, List[WernerStateEntanglement]] = {} #初始化资源情况
            while next_hop!=dest and next_hop!=None: #开始纠缠分发
                j1 = i1+1
                # print(len(path))
                # print("i=",i1)
                # print("j=",j1)

                current = path[i1]
                next_hop = path[j1]
                l1 = len(current.qroute_table[next_hop]) #找纠缠链路
                if l1>=m:   #资源充足，不用恢复路径
                    ql: List[WernerStateEntanglement]=current.qroute_table.get(next_hop,None)
                    quantum_source[next_hop] = ql  # 找到有m份资源的候选节点
                    for item, item2 in quantum_source.items():
                        for n1 in range(0, len(quantum_source[item])):
                            for n2 in range(n1 + 1, len(quantum_source[item])):
                                if quantum_source[item][n1].w < quantum_source[item][n2].w:
                                    quantum_source[item][n1], quantum_source[item][n2] = quantum_source[item][n2], \
                                                                                       quantum_source[item][n1]
                    # print(quantum_source)  ##得到最终的资源结果,对每个节点的纠缠资源按照保真度降序排列
                    if m == 1:
                        if current == src:
                            epr = quantum_source[next_hop][0]
                            next_hop.qroute_table[current].remove(epr)
                            current.qroute_table[next_hop].remove(epr)  #选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表

                            current.memories[0].read(epr)
                            next_hop.memories[0].read(epr)
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = epr
                        else:
                            z = quantum_source[next_hop][0]
                            current.memories[0].read(z)
                            next_hop.memories[0].read(z)
                            epr = epr.swapping(z)
                            del current.qroute_table[next_hop][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                            del next_hop.qroute_table[current][0]
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = []  # 将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                            if epr == -1:
                                epr = None
                                    # print("fail to swap")
                                return  ##纠缠交换失败，资源作废，请求失败
                        # i=i+1
                        # continue
                            # print("epr is",epr)
                    elif m > 1:

                        while (len(current.qroute_table[next_hop]) >= m):
                            m1 = m - 1
                            for i in range(0, m1):  # range不包括结束的节点
                                current.memories[0].read(quantum_source[next_hop][m - i - 1])
                                next_hop.memories[0].read(quantum_source[next_hop][m - i - 1])
                                # s1 = quantum_source[next_hop][i]

                                # s = quantum_source[next_hop][num+i]
                                current.memories[0].read(quantum_source[next_hop][m - i - 2])
                                next_hop.memories[0].read(
                                    quantum_source[next_hop][m - i - 2])  ##将最终选定的节点中需要用到的资源从节点取出，并在qroutetable上删除
                                # current.qroute_table[next_hop].remove(quantum_source[next_hop][m-i-1])
                                # next_hop.qroute_table[current].remove(quantum_source[next_hop][m-i-1])
                                # current.qroute_table[next_hop].remove(quantum_source[next_hop][num+i])
                                # next_hop.qroute_table[current].remove(quantum_source[next_hop][num+i])
                                a = quantum_source[next_hop][m - i - 2]
                                b = quantum_source[next_hop][m - i - 1]

                                quantum_source[next_hop][m - i - 2] = a.distillation(b)
                                del next_hop.qroute_table[current][m - i - 1]
                                del current.qroute_table[next_hop][m - i - 1]
                                epr1 = quantum_source[next_hop][m - i - 2]
                                if epr1 == -1:
                                    del next_hop.qroute_table[current][m - i - 2]
                                    del current.qroute_table[next_hop][m - i - 2]
                                    # print("fail to distill")
                                    break
                            # 此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤
                            if epr1 == -1:
                                epr1 = None
                                if (len(current.qroute_table[next_hop]) < m):
                                    return
                                continue #纯化失败，重新找资源进行纯化
                                # current.memories[0].write(quantum_source[next_hop][0])
                                # next_hop.memories[0].write(quantum_source[next_hop][0])
                            if current == src:
                                epr = epr1
                                del current.qroute_table[next_hop][0]
                                del next_hop.qroute_table[current][0]
                                for key in lock_table:
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                            else:
                                # print(epr)
                                epr = epr.swapping(epr1)  # epr 指向延长的纠缠资源
                                del current.qroute_table[next_hop][0]
                                del next_hop.qroute_table[current][0]
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                    # print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                                    return
                            break #如果成功建立一个纠缠资源，就不再需要在路径上找资源了，跳出循环
                    i1 = i1 + 1  # 延长了一跳，去下一跳
                    if epr.fidelity < Fth:
                        epr = -1
                        m = m + m
                        break

                elif l1<m: #资源不足，找恢复路径
                    sd = (current, next_hop)
                    l3 = math.inf
                    rpath4 = []
                    pathlist = RecoveryPath.get(sd, None)
                    if pathlist == None:
                        return  # 无恢复路径
                    for rpath3 in pathlist:
                        if len(rpath3) < l3 and len(current.qroute_table[rpath3[1]]) >= m:
                            l3 = len(rpath3)
                            rpath4 = rpath3.copy()
                            # 找有足够资源且最短的一条作为恢复路径
                    if rpath4 == []:
                        return  # 找不到恢复路径，且主路径没资源，请求失败
                    for node in rpath4:
                        if node == current or node == next_hop:
                            rpath4.remove(node)
                    rpath4.reverse()
                    for node in rpath4:
                        path.insert(j1, node)  # 将恢复路径插入原路径
                    break
            if epr==-1:
                epr = None
                continue
            while(next_hop!=dest and next_hop!=None): #接着用恢复路径完成请求
                j1=i1+1
                current = path[i1]
                next_hop = path[j1]
                l2 = len(current.qroute_table[next_hop])
                if l2>=m:
                    ql:List[WernerStateEntanglement]=current.qroute_table.get(next_hop,None)
                    quantum_source[next_hop] = ql  # 找到有m份资源的候选节点
                    for item, item2 in quantum_source.items():
                        for n1 in range(0, len(quantum_source[item])):
                            for n2 in range(n1 + 1, len(quantum_source[item])):
                                if quantum_source[item][n1].w < quantum_source[item][n2].w:
                                    quantum_source[item][n1], quantum_source[item][n2] = quantum_source[item][n2], \
                                                                                           quantum_source[item][n1]
                        # print(quantum_source)  ##得到最终的资源结果,对每个节点的纠缠资源按照保真度降序排列
                    if m == 1:
                        if current == src:
                            epr = quantum_source[next_hop][0]
                            next_hop.qroute_table[current].remove(epr)
                            current.qroute_table[next_hop].remove(epr)  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                            current.memories[0].read(epr)
                            next_hop.memories[0].read(epr)
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = epr
                        else:
                            z = quantum_source[next_hop][0]
                            current.memories[0].read(z)
                            next_hop.memories[0].read(z)
                            epr = epr.swapping(z)
                            del current.qroute_table[next_hop][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                            del next_hop.qroute_table[current][0]
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  # 将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                            if epr == -1:
                                epr = None
                                        # print("fail to swap")
                                return  ##纠缠交换失败，资源作废，请求失败
                                # print("epr is",epr)
                    elif m > 1:
                        while (len(current.qroute_table[next_hop]) >= m):
                            m1 = m - 1
                            for i in range(0, m1):  # range不包括结束的节点
                                current.memories[0].read(quantum_source[next_hop][m - i - 1])
                                next_hop.memories[0].read(quantum_source[next_hop][m - i - 1])
                                    # s1 = quantum_source[next_hop][i]

                                    # s = quantum_source[next_hop][num+i]
                                current.memories[0].read(quantum_source[next_hop][m - i - 2])
                                next_hop.memories[0].read(
                                        quantum_source[next_hop][m - i - 2])  ##将最终选定的节点中需要用到的资源从节点取出，并在qroutetable上删除
                                    # current.qroute_table[next_hop].remove(quantum_source[next_hop][m-i-1])
                                    # next_hop.qroute_table[current].remove(quantum_source[next_hop][m-i-1])
                                    # current.qroute_table[next_hop].remove(quantum_source[next_hop][num+i])
                                    # next_hop.qroute_table[current].remove(quantum_source[next_hop][num+i])
                                a = quantum_source[next_hop][m - i - 2]
                                b = quantum_source[next_hop][m - i - 1]

                                quantum_source[next_hop][m - i - 2] = a.distillation(b)
                                del next_hop.qroute_table[current][m - i - 1]
                                del current.qroute_table[next_hop][m - i - 1]
                                epr1 = quantum_source[next_hop][m - i - 2]
                                if epr1 == -1:
                                    del next_hop.qroute_table[current][m - i - 2]
                                    del current.qroute_table[next_hop][m - i - 2]
                                        # print("fail to distill")
                                    break


                                # 此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤
                            if epr1 == -1:
                                epr1 = None
                                if (len(current.qroute_table[next_hop]) < m):  # 资源不足，无法满足请求
                                    return
                                continue
                                    # current.memories[0].write(quantum_source[next_hop][0])
                                    # next_hop.memories[0].write(quantum_source[next_hop][0])
                            if current == src:
                                epr = epr1
                                del current.qroute_table[next_hop][0]
                                del next_hop.qroute_table[current][0]
                                for key in lock_table:
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = epr
                            else:
                                    # print(epr)
                                epr = epr.swapping(epr1)  # epr 指向延长的纠缠资源
                                del current.qroute_table[next_hop][0]
                                del next_hop.qroute_table[current][0]
                                for key, value in lock_table.items():
                                    if (key.src == src and key.dest == dest):
                                        lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                                if epr == -1:
                                        # print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                                    return
                            break #如果成功建立一个纠缠资源，就不再需要在路径上找资源了，跳出循环
                else:
                    return #候选路径无资源

                i1=i1+1 #延长了一跳，去下一跳
                if epr.fidelity < Fth:
                    m=m+m
                    break
            if(epr.fidelity>Fth):
                print(epr.fidelity)
                return True
        if m>=16:
            return #需要的资源过多，无法实现纯化



class QCASTRouteAlgorithm2(RouteImpl):
    """
    This is the dijkstra route algorithm implement
    """

    INF = math.inf

    def __init__(self, name: str = "qcast",
                 metric_func: Callable[[Union[QuantumChannel, ClassicChannel]], float] = None) -> None:
        """
        Args:
            name: the routing algorithm's name
            metric_func: the function that returns the metric for each channel.
                The default is the const function m(l)=1
        """
        self.name = name
        self.route_table = {}
        if metric_func is None:
            self.metric_func = lambda _: 1
        else:
            self.metric_func = metric_func

    def build(self, nodes: List[QNode], channels: List[Union[QuantumChannel, ClassicChannel]]):

        for n in nodes:
            selected = []
            unselected = [u for u in nodes]

            d = {}
            for nn in nodes:
                if nn == n:
                    d[n] = [0, []]
                else:
                    d[nn] = [self.INF, [nn]]

            while len(unselected) != 0:
                ms = unselected[0]
                mi = d[ms][0]

                for s in unselected:
                    if d[s][0] < mi:
                        ms = s
                        mi = d[s][0]

                # d[ms] = [d[ms][0], d[ms][1]]
                selected.append(ms)
                unselected.remove(ms)

                for link in channels:
                    if ms not in link.node_list:
                        continue
                    if len(link.node_list) < 2:
                        raise NetworkRouteError("broken link")
                    idx = link.node_list.index(ms)
                    idx_s = 1 - idx
                    s = link.node_list[idx_s]
                    if s in unselected and d[s][0] > d[ms][0] + self.metric_func(link):
                        d[s] = [d[ms][0] + self.metric_func(link), [ms] + d[ms][1]]

            for nn in nodes:
                d[nn][1] = [nn] + d[nn][1]
            self.route_table[n] = d

    def minlen(self,current:QNode,dest:QNode):
        """
        find the length of shortest path to build the E2E entanglement

        """

        n1 = self.route_table.get(current,None)
        # n1: Dict[QNode, List[float, List[QNode]]] = self.route_table.get(current, None)
        for key in n1:
            if key == dest:
                minlen = n1[key][0]
        return minlen

    def query(self, src: QNode, dest: QNode,lock_table: Dict[Request,WernerStateEntanglement],Fth:float) -> List[Tuple[float, QNode, List[QNode]]]:
        """
        query the metric, nexthop and the path

        Args:
            src: the source node
            dest: the destination node

        Returns:
            A list of route paths. The result should be sortted by the priority.
            The element is a tuple containing: metric, the next-hop and the whole path.
        """
        current = src
        fidelity=0.97
        # find the candidate of next hop
        RecoveryPath = {} #sd元组对key的恢复路径集合
        # print("qroute_table is",current.qroute_table)
        epr: WernerStateEntanglement = None  # 标记延长的纠缠资源 ,初始化epr
        m = 1

            # print("m =",m) ##得到m，m是在寻找下一跳时需要的资源
        ls: Dict[QNode, List[float, List[QNode]]] = self.route_table.get(src, None)
        if ls is None:
            return
        le = ls.get(dest, None)
        if le is None:
            return
        path: List[QNode] = le[1]
        path = path.copy()
        path.reverse() #主路径path，接下来找恢复路径
        path1 = path.copy()

        l = len(path)-1
        for i in range(0,l):
            if i+1<=l:
                for j in range(i+1,l+1):
                    nei1 = set()
                    nei2 = set()
                    current = path[i]
                    dest1 = path[j]
                    n1 = self.route_table.get(current,None)
                    n2 = self.route_table.get(dest1,None)
                    for key in n2:
                        nei2.add(key)
                    for key in n1:
                        if (n1[key][0] == 1):
                            nei1.add(key)
                            # find the neighbors of current
                    for key in n1:
                            # print("error")
                        if (key == dest1):
                            minlen = n1[key][0]  # 源节点和目的节点之间的距离
                    C = nei1 & nei2  # current的邻居中能到dest1的节点集合
                    for z in path:
                        C.discard(z) #把路径上的节点去掉
                    sd: Tuple(QNode, QNode) = (current, dest1)
                    Recovery = []
                    for k in C:
                        pathid=0
                        rpath = n2.get(k,None) #点k到点dest1的最短路径
                        for node in rpath[1]:
                            if node==current:
                                pathid=1
                                break
                        if pathid==1:
                            pathid=0
                            continue
                        if rpath[0]<3:
                            rpath1:List[QNode] = rpath[1]
                            rpath1 = rpath1.copy()
                            rpath1.insert(0,current)
                            Recovery.append(rpath1)

                    RecoveryPath[sd]=Recovery
                        #找到恢复路径，开始建立纠缠链接,先生成资源，再进行纠缠交换
        while (m<16):  #当生成纠缠链接保真度不足时会回退到这里

            #if (epr == -1): #需要更多的纠缠资源
            entanglementlist: List[WernerStateEntanglement]=[]
            quantum_source: Dict[QNode, List[WernerStateEntanglement]] = {}
            path=path1.copy()
            l = len(path)-1
            i1=0
            current=path[0]
            next_hop=path[0]

            while next_hop!=dest and next_hop!=None: #开始纠缠分发，

                j1 = i1+1
                current = path[i1]
                next_hop = path[j1]
                current.qroute_table[next_hop]=[]
                next_hop.qroute_table[current]=[]
                for k in range(0, m*2):#生成纠缠资源，m大于1时，需要将资源进行纠缠纯化
                    name_id=1
                    k1 = get_rand(0, 1)
                    if k1 < 0.5:    # 随机生成纠缠链接，概率0.5
                        name_id = k + 1
                        e1 = WernerStateEntanglement(name=f"e{name_id}", fidelity=fidelity)  # 生成m份纠缠资源
                        current.memories[0].write(e1) #写入
                        next_hop.memories[0].write(e1)
                        current.qroute_table[next_hop].append(e1)                   #修改qroutetable，将资源放进去
                        next_hop.qroute_table[current].append(e1)
                if len(current.qroute_table[next_hop]) < m: #资源生成数量不足，需要恢复路径

                    sd=(current,next_hop)
                    l3 = math.inf
                    rpath4 = []
                    pathlist= RecoveryPath.get(sd, None) #得到恢复路径队列
                    if pathlist==None:
                        return                           #无恢复路径
                    for rpath3 in pathlist:
                        if len(rpath3) < l3:
                            l3 = len(rpath3)
                            rpath4 = rpath3.copy()    #选择其中最短的恢复路径
                        # if rpath4 == []:
                        #     return  # 找不到恢复路径，且主路径没资源，请求失败
                    for node in rpath4:
                        if node == current or node == next_hop:
                            rpath4.remove(node)
                    rpath4.reverse()
                    for node in rpath4:
                        path.insert(j1, node)  # 将恢复路径插入原路径

                    current = path[i1]
                    next_hop = path[j1]
                    current.qroute_table[next_hop]=[]
                    next_hop.qroute_table[current]=[]
                    for k in range(0, m * 2):  # 生成纠缠资源，m大于1时，需要将资源进行纠缠纯化
                        name_id = 1
                        k1 = get_rand(0, 1)
                        if k1 < 0.5:  # 随机生成纠缠链接，概率0.5
                            name_id = k + 1
                            e1 = WernerStateEntanglement(name=f"e{name_id}", fidelity=fidelity)  # 生成m份纠缠资源
                            current.memories[0].write(e1)  # 写入
                            next_hop.memories[0].write(e1)
                            current.qroute_table[next_hop].append(e1)  # 修改qroutetable，将资源放进去
                            next_hop.qroute_table[current].append(e1)
                        # else:
                        #     return #恢复路径上生成资源也失败，视为纠缠生成失败,但是会继续进行纠缠资源的生成
                    if len(current.qroute_table[next_hop])<m:#检查资源生成数量是否达到m个
                        return
                    elif len(current.qroute_table[next_hop]) >= m:
                        ql: List[WernerStateEntanglement] = current.qroute_table.get(next_hop, None)
                        quantum_source[next_hop] = ql
                elif len(current.qroute_table[next_hop])>=m:
                    ql: List[WernerStateEntanglement] = current.qroute_table.get(next_hop, None)
                    quantum_source[next_hop] = ql
                    # 找到有m份资源的候选节点


# 资源生成结束，开始纠缠交换

                if m == 1:
                    if current == src:
                        epr = quantum_source[next_hop][0]
                        next_hop.qroute_table[current].remove(epr)
                        current.qroute_table[next_hop].remove(epr)  #选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表

                        current.memories[0].read(epr)
                        next_hop.memories[0].read(epr)
                        for key, value in lock_table.items():
                            if (key.src == src and key.dest == dest):
                                lock_table[key] = epr
                    else:
                        z = quantum_source[next_hop][0]
                        current.memories[0].read(z)
                        next_hop.memories[0].read(z)
                        epr = epr.swapping(z)
                        del current.qroute_table[next_hop][0]  ##选定的下一跳资源将会从量子资源表和节点中删除，并写入锁定资源表
                        del next_hop.qroute_table[current][0]
                        for key, value in lock_table.items():
                            if (key.src == src and key.dest == dest):
                                lock_table[key] = []  # 将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                        if epr == -1:
                            epr = None                                    # print("fail to swap")
                            return  ##纠缠交换失败，资源作废，请求失败
                            # print("epr is",epr)
                elif m > 1:
                    while (len(current.qroute_table[next_hop]) >= m):
                        m1 = m - 1
                        for i in range(0, m1):  # range不包括结束的节点
                            current.memories[0].read(quantum_source[next_hop][m - i - 1])
                            next_hop.memories[0].read(quantum_source[next_hop][m - i - 1])
                            # s1 = quantum_source[next_hop][i]

                            # s = quantum_source[next_hop][num+i]
                            current.memories[0].read(quantum_source[next_hop][m - i - 2])
                            next_hop.memories[0].read(
                                quantum_source[next_hop][m - i - 2])  ##将最终选定的节点中需要用到的资源从节点取出，并在qroutetable上删除
                            # current.qroute_table[next_hop].remove(quantum_source[next_hop][m-i-1])
                            # next_hop.qroute_table[current].remove(quantum_source[next_hop][m-i-1])
                            # current.qroute_table[next_hop].remove(quantum_source[next_hop][num+i])
                            # next_hop.qroute_table[current].remove(quantum_source[next_hop][num+i])
                            a = quantum_source[next_hop][m - i - 2]
                            b = quantum_source[next_hop][m - i - 1]

                            quantum_source[next_hop][m - i - 2] = a.distillation(b)
                            del next_hop.qroute_table[current][m - i - 1]
                            del current.qroute_table[next_hop][m - i - 1]
                            epr1 = quantum_source[next_hop][m - i - 2]
                            if epr1 == -1:
                                del next_hop.qroute_table[current][m - i - 2]
                                del current.qroute_table[next_hop][m - i - 2]
                                # print("fail to distill")
                                break

                        # 此时，当纠缠纯化失败后，我没有将所有的资源都清除，而是只清除了纯化失败的那两个，纠缠纯化失败退回到寻找下一跳的步骤
                        if epr1 == -1:
                            epr1 = None
                            if (len(current.qroute_table[next_hop]) < m):  # 资源不足，无法满足请求
                                return
                            continue
                            # current.memories[0].write(quantum_source[next_hop][0])
                            # next_hop.memories[0].write(quantum_source[next_hop][0])
                        if current == src:
                            epr = epr1
                            del current.qroute_table[next_hop][0]
                            del next_hop.qroute_table[current][0]
                            for key in lock_table:
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = epr
                        else:
                            # print(epr)
                            epr = epr.swapping(epr1)  # epr 指向延长的纠缠资源
                            del current.qroute_table[next_hop][0]
                            del next_hop.qroute_table[current][0]
                            for key, value in lock_table.items():
                                if (key.src == src and key.dest == dest):
                                    lock_table[key] = []  ##将保留的第一跳纠缠资源解锁，但是这里在后续每一次都删除了一遍，理论上不需要
                            if epr == -1:
                                # print("fail to swap")  ##纠缠交换失败，应该跳出循环,本次寻找路径失败
                                return
                i1=i1+1

                if epr.fidelity < Fth:
                    m=m+m
                    break
            if(epr.fidelity>Fth):
                return True
        if m>=16:
            return #需要的资源过多，无法实现纯化
