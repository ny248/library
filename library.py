class queue_elem:
    def __init__(self):
        self.next = None
        self.prev = None
        self.value = None
class queue:
    def __init__(self):
        self.begin = queue_elem()
        self.end = queue_elem()
        self.begin.next = self.end
        self.end.pref = self.begin
    def empty(self):
        return self.begin.next == self.end
    def add_left(self, value):
        elem = queue_elem()
        elem.prev = self.begin
        elem.next = self.begin.next
        self.begin.next.prev = elem
        self.begin.next = elem
        elem.value = value
    def add_right(self, elem):
        elem = queue_elem()
        elem.prev = self.end.prev
        elem.next = self.end
        self.end.prev.next = elem
        self.end.prev = elem
        elem.value = value
    def pop_left(self):
        if self.empty():
            return None
        ret = self.begin.next.value
        self.begin.next = self.begin.next.next
        del self.begin.next.prev
        self.begin.next.prev = self.begin
        return ret
    def pop_right(self):
        if self.empty():
            return None
        ret = self.end.prev.value
        self.end.prev = self.end.prev.prev
        del self.end.prev.next
        self.end.prev.next = self.end
        return ret

class mod_int: #素数を法とした自然数, または有理数を扱うことのできるクラス
    def __init__(self, value, mod):
        self.value = value % mod
        self.mod = mod
    def __neg__(self):
        return mod_int(-self.value, self.mod)
    def __add__(self, other):
        if type(other) == mod_int:
            return mod_int(self.value + other.value, self.mod)
        elif type(other) == int:
            return mod_int(self.value + other, self.mod)
        raise TypeError()
    def __sub__(self, other):
        return self + (-other)
    def __mul__(self, other):
        if type(other) == mod_int:
            return mod_int(self.value * other.value, self.mod)
        elif type(other) == int:
            return mod_int(self.value * other, self.mod)
        raise TypeError()
    def __truediv__(self, other):
        if type(other) in [mod_int, int]:
            return self * other ** (self.mod - 2)
        raise TypeError()
    def __pow__(self, other):
        if type(other) != int:
            raise TypeError()
        cur, ret = self.value, 1
        while other > 0:
			if other % 2:
                ret = (ret * cur) % self.mod
            other //= 2
            cur = (cur ** 2) % self.mod
        return mod_int(ret, self.mod)
    def __repr__(self):
        return str(self.value)
A = mod_int(3, 998244353)
print(A / 6)

# 入力:大きさ2 ** N の配列 A[i]
# 出力:大きさ2 ** N の配列 B[i] B[i] は i | j = i　となるような全てのjに対してのA[j]の和
def subset_zeta(A, N):
    B = A.copy()
    mask = 1
    for i in range(N):
        for j in range(2 ** N):
            if j & mask == 0:
                B[j | mask] += B[j]
        mask *= 2
    return B
subset_zeta([2,5,3,7],2)
#逆変換
def subset_mobius(A, N):
    B = A.copy()
    mask = 1
    for i in range(N):
        for j in range(2 ** N):
            if j & mask:
                B[j] -= B[j ^ mask]
        mask *= 2
    return B
print(subset_mobius([2,7,5,17],2))

# 入力:大きさ2 ** N の配列 A[i]
# 出力:大きさ2 ** N の配列 B[i] B[i] は A[k] * exp(-2 * j * math.pi * i * k / N)のjについての和(j は虚数単位)
# 参考 : https://ja.wikipedia.org/wiki/%E9%AB%98%E9%80%9F%E3%83%95%E3%83%BC%E3%83%AA%E3%82%A8%E5%A4%89%E6%8F%9B
import math, cmath
def fourier(A, N):
    if N == 0:
        return A
    P = 1 << (N - 1)
    ret = [0 for i in range(2 * P)]
    for r in range(2):
        f = [0 for i in range(P)]
        for p in range(P):
            for q in range(2):
                f[p] += A[q * P + p] * cmath.exp(-2 * math.pi * complex(0, 1) * r * q / 2)
            f[p] *= cmath.exp(-math.pi * complex(0, 1) * p * r / P)
        f = fourier(f, N - 1)
        for p in range(len(f)):
            ret[2 * p + r] = f[p]
    return ret
def inv_fourier(A, N):
    coeff = 1 / (1 << N)
    B = [i.conjugate() for i in A]
    B = fourier(B, N)
    for i in range(1 << N):
        B[i] = B[i].conjugate() * coeff
    return B
inv_fourier(fourier([1,2,3,4,5,6,7,8],3),3)

class segment_tree:
    def __init__(self,A,unit,func):
        self.N=1
        self.unit=unit
        self.func=func
        while self.N<len(A):
            self.N*=2
        self.memo=[unit for i in range(self.N*2-1)]
        for i in range(len(A)):
            self.memo[i+self.N-1]=A[i]
        for i in range(self.N-2,-1,-1):
            self.memo[i]=func(self.memo[2*i+1],self.memo[2*i+2])
    def change(self,index,value):
        i=index+self.N-1
        self.memo[i]=value
        while i>0:
            i=(i-1)//2
            self.memo[i]=self.func(self.memo[2*i+1],self.memo[2*i+2])
    def get(self,begin,end):
        b,e=begin+self.N-1,end+self.N-1
        ans=self.unit
        while b<e:
            if b%2==0:
                ans=self.func(ans,self.memo[b])
            b=b//2
            if e%2==0:
                ans=self.func(ans,self.memo[e-1])
            e=(e-1)//2
        return ans
A=segment_tree([1,6,1,4,6,3,3],0,lambda a,b:a+b)
print(A.get(0,6))
A.change(2,9)
print(A.get(1,4))



def kadane(A):
    N=len(A)
    S=[0 for i in range(N)]
    S,ans=A[0],A[0]
    for i in range(1,N):
        S=max(S,0)+A[i]
        ans=max(S,ans)
    return ans
kadane([-5,-1,6,4,9,-6,-7])

def divisors(n):
    ans=[]
    for i in range(1,n+1):
        if n%i==0:
            ans.append(i)
        if i**2>n:
            break
    k=len(ans)
    for i in range(k-1,-1,-1):
        if n//ans[i]>ans[k-1]:
            ans.append(n//ans[i])
    return ans
print(divisors(5))

def gcd(a,b):
    return gcd(a%b,b%a) if min(a,b,abs(a-b)) else max(a,b)
def lcm(a,b):
    return a*b//gcd(a,b)

class unionFind:
    def __init__(self,N):
        self.N=N
        self.parent=[-1 for i in range(N)]
        self.size=[1 for i in range(N)]
    def find(self,x):
        path=[x]
        while self.parent[path[-1]]!=-1:
            path.append(self.parent[path[-1]])
        for i in path[:-1]:
            self.parent[i]=path[-1]
        return path[-1]
    def unite(self,x,y):
        roots=sorted([self.find(x),self.find(y)],key=lambda _:self.parent[_])
        if roots[0]!=roots[1]:
            self.parent[roots[0]]=roots[1]
            self.size[roots[1]]+=self.size[roots[0]]

U=unionFind(4)
U.unite(1,2)
U.unite(1,3)
U.unite(1,3)
for i in range(4):
    print(U.find(i))

class queue:
    def __init__(self,A):
        self.free=[]
        self.data=A
        self.next=[i+1 for i in range(len(A)-1)]+[-1]*(len(A)>0)
        self.head=0
        self.tail=len(self.data)-1
        self.N=len(A)
    def push(self,value):
        if len(self.free):
            ind=self.free[-1]
            self.free.pop()
        else:
            ind=len(self.data)
            self.data.append(0)
            self.next.append(0)
        self.data[ind]=value
        self.next[ind]=-1
        if self.N>0:
            self.next[self.tail]=ind
        else:
            self.head=ind
        self.tail=ind
        self.N+=1
    def front(self):
        if self.N:
            return self.data[self.head]
        else:
            raise
    def pop(self):
        if self.N:
            ret=self.front()
            self.free.append(self.head)
            self.head=self.next[self.head]
            self.N-=1
            return ret
        else:
            raise
q=queue([])

class priorityQueue:
    def __init__(self,key,A=[]):
        self.data=[]
        self.N=0
        self.key=key
        for i in A:
            self.push(i)
    def children(self,index):
        return [i for i in range(index*2+1,min(self.N,index*2+3))]
    def parent(self,index):
        return (index-1)//2
    def push(self,value):
        self.N+=1
        self.data.append(value)
        cur=self.N-1
        while cur>0:
            if self.key(self.data[(cur-1)//2],self.data[cur]):
                self.data[cur],self.data[(cur-1)//2]=self.data[(cur-1)//2],self.data[cur]
            cur=(cur-1)//2
    def top(self):
        return self.data[0]
    def pop(self):
        ret=self.data[0]
        self.data[0]=self.data[self.N-1]
        self.data.pop()
        self.N-=1
        cur=0
        while len(self.children(cur)):
            c=self.children(cur)
            if all([self.key(self.data[i],self.data[cur]) for i in c]):
                break
            target=(len(c)==2 and self.key(self.data[c[0]],self.data[c[1]]))
            self.data[cur],self.data[c[target]]=self.data[c[target]],self.data[cur]
            cur=c[target]
        return ret

q=priorityQueue(lambda a,b:a>=b,[1,6,1,3,6,2,24,52,1,1252,521])
for i in range(11):
    print(q.pop())
    q.push(4)
    print(q.data)

import heapq
class node:
    def __init__(self, con):
        self.con = con#つながっている辺の行き先と重みの配列
class directedWeightedGraph:
    def __init__(self, N, edges = []):#edges:[from,to,weight]の配列
        self.N = N
        self.nodes = [node([]) for i in range(N)]
        for i in edges:
            self.nodes[i[0]].append(i[1:])
    def connect(self, f, t, weight = 0):
        self.nodes[f].con.append([t, weight])
    def dijkstra(self, start):#startから各点への距離を返す
        ans = [float('inf') for i in range(self.N)]
        watch = []
        heapq.heappush(watch, [0, start])
        while watch:
            cur = heapq.heappop(watch)#今見ているノード
            if ans[cur[1]] >= cur[0]:
                ans[cur[1]] = cur[0]
                for i in self.nodes[cur[1]].con:
                    heapq.heappush(watch, [i[1] + ans[cur[1]], i[0]])
        return ans
G=directedWeightedGraph(3)
G.connect(0,1,3)
G.connect(1,2,4)
G.connect(0,2,19)
print(G.dijkstra(0))

def p(a,p,mod):#素数体上の累乗
    rem=p%(mod-1)
    temp=a
    ans=1
    while rem:
        if rem%2:
            ans=(ans*temp)%mod
        temp=(temp**2)%mod
        rem//=2
    return ans
p(3,6,5)

import heapq
heap = []
data = [[1, 'J'], [4, 'N'], [3, 'H'], [2, 'O']]
for item in data:
    heapq.heappush(heap, item)
while heap:
    print(heapq.heappop(heap)[1])

def knapsack(items,capacity):#itemsは[重さ、価値]の配列
    N=len(items)
    dp=[[0 for i in range(capacity+1)] for j in range(N+1)]
    for i in range(N):
        for j in range(min(capacity+1,items[i][0])):
            dp[i+1][j]=dp[i][j]
        for j in range(items[i][0],capacity+1):
            dp[i+1][j]=max(dp[i][j],dp[i][j-items[i][0]]+items[i][1])
    return dp[N][capacity]

def euclid(a,b,c):#ax+by=cの解を一つ返す 解が存在しない時は-1をかえす
    print(a,b,c)
    if a==b==0:
        return [[1,1],-1][c!=0]
    if a<0 or b<0:
        temp=euclid(abs(a),abs(b),c)
        return -1 if temp==-1 else [[1,-1][a<0]*temp[0],[1,-1][b<0]*temp[1]]
    if a==0:
        return -1 if c%b else [0,c//b]
    if b==0:
        return -1 if c%a else [c//a,0]
    if b>a:
        temp=euclid(a,b%a,c)
        print(temp)
        return -1 if temp==-1 else [temp[0]-(b//a)*temp[1],temp[1]]
    else:
        temp=euclid(a%b,b,c)
        print(temp)
        return -1 if temp==-1 else [temp[0],temp[1]-(a//b)*temp[0]]
print(euclid(13,14,4))

def insert(A):
    N=len(A)
    ans=0
    for i in range(1,N):
        j=i
        while j>0 and A[j]<A[j-1]:
            A[j],A[j-1]=A[j-1],A[j]
            j-=1
            ans+=1
    return ans
def merge(A,l,r):
    if r-l<2:
        return 0
    m=(r+l)//2
    ans=merge(A,l,m)+merge(A,m,r)
    B=[]
    pt=[l,m]
    while pt[0]<m or pt[1]<r:
        if pt[0]==m or (pt[1]<r and A[pt[0]]>A[pt[1]]):
            B.append(A[pt[1]])
            pt[1]+=1
            ans+=m-pt[0]
        else:
            B.append(A[pt[0]])
            pt[0]+=1
    A[l:r]=B
    return ans
print(merge([8,3,6,5,2,4,1,9,7]))

def factorization(a):
    i=2
    ans=[]
    while i*i<=a:
        while a%i==0:
            ans.append(i)
            a//=i
        i+=1
    if a>1:
        ans.append(a)
    return sorted(ans)
print(factorization(20385))

def next_combination(A):#0と1からなる配列でcombinationを表す 辞書順で次の配列を返す 終わったら-1を返す O(N)
    N=len(A)
    K=sum(A)
    if N==1 or K==0 or K==N:
        return -1
    B=next_combination(A[1:])
    if B==-1:
        if A[0]==0:
            return [1]+[0]*(N-K)+[1]*(K-1)
        else:
            return -1
    else:
        return [A[0]]+B
P=[0,0,0,1,1]
while P!=-1:
    P=next_combination(P)
    print(P)
