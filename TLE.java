import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.abs;
import static java.lang.Math.ceilDiv;
import java.io.*;
import java.util.*;

public class TLE {
    public static PrintWriter out = new PrintWriter(new BufferedOutputStream(System.out));
    static FastReader in = new FastReader();
    // static final int MOD = (int) 1e9 + 7;
    static final int MOD = 998244353;
    static List<List<Integer>> adj;
    static final int MAX_LOG = 20;
    static final int N = (int) 2e5 + 1;
    static int[][] par = new int[N][MAX_LOG];
    static int[] depth = new int[N];
    static long maxCost = 0;

    /*
     * @Sulabh Ambule
     */

    public static void main(String[] Hi) throws IOException {
        int T = 1;
        while (T-- > 0) {
            Code_Subh();
        }
        out.close();
    }

    static void Code_Subh() throws IOException {
        int n = in.nextInt();
        int[] a = inputIntArr(n);
        adj = new ArrayList<>();
        for (int i = 0; i <= n; i++) {
            adj.add(new ArrayList<>());
        }

        for (int i = 0; i < n - 1; i++) {
            int u = in.nextInt();
            int v = in.nextInt();
            adj.get(u).add(v);
            adj.get(v).add(u);
        }
        long sum = 0;
        for (int i : a) {
            sum += i;
        }
        int[] subTreeSum = new int[n + 1];
        int[] ans1 = { 0 };
        dfs(1, 0, 0, ans1, a, subTreeSum);
        maxCost = ans1[0];
        dfs3(1, 0, maxCost, subTreeSum, sum);
        out.println(maxCost);
    }

    static void dfs(int node, int parent, int d, int[] ans1, int[] a, int[] subTreeSum) {
        subTreeSum[node] += a[node - 1];
        ans1[0] += d * a[node - 1];
        for (int child : adj.get(node)) {
            if (parent == child)
                continue;
            dfs(child, node, d + 1, ans1, a, subTreeSum);
            subTreeSum[node] += subTreeSum[child];
        }
    }

    static void dfs3(int node, int parent, long currCost, int[] subTreeSum, long total) {
        maxCost = max(maxCost, currCost);
        for (int adjNode : adj.get(node)) {
            if (adjNode != parent) {
                // re Roting the tree here
                long newCost = currCost - subTreeSum[adjNode] + (total - subTreeSum[adjNode]);
                dfs3(adjNode, node, newCost, subTreeSum, total);
            }
        }

    }

    // static void dfs(int node, int parent, int[] subTreeSum) {
    // depth[node] = 1 + depth[parent];
    // par[node][0] = parent;
    // for (int j = 0; j < MAX_LOG; j++) {
    // par[node][j] = par[par[node][j - 1]][j - 1];
    // }
    // for (int adjNode : adj.get(node)) {
    // if (adjNode != parent) {
    // dfs(adjNode, node, subTreeSum);
    // subTreeSum[node] += subTreeSum[adjNode];
    // }
    // }
    // }

    // static int LCA(int u, int v) {
    // if (u == v)
    // return u;
    // if (depth[u] < depth[v]) {
    // int t = u;
    // u = v;
    // v = t;
    // }
    // int d = depth[u] - depth[v];
    // for (int i = MAX_LOG - 1; i >= 0; i--) {
    // if (((1 << i) & d) != 0) {
    // u = par[u][i];
    // }
    // }
    // for (int j = MAX_LOG - 1; j >= 0; j--) {
    // if (par[u][j] != par[v][j]) {
    // u = par[u][j];
    // v = par[v][j];
    // }
    // }
    // return (u != v) ? par[u][0] : u;
    // }

    /*------------------------------------------------------------------------------------------------- */

    static class Pair implements Comparable<Pair> {
        long first;
        int second;

        Pair(long f, int s) {
            this.first = f;
            this.second = s;
        }

        @Override
        public int compareTo(Pair other) {
            return (int) (this.first - other.first);
        }
    }

    static int gcd(int a, int b) {
        if (a == 0)
            return b;
        return gcd(b % a, a);
    }

    static int lcm(int a, int b) {
        return Math.abs(a * b) / gcd(a, b);
    }

    static boolean isPrime(long arr) {
        if (arr <= 1)
            return false;
        for (int i = 2; i <= Math.sqrt(arr); i++) {
            if (arr % i == 0)
                return false;
        }
        return true;
    }

    static void reverse(int[] a, int l, int r) {
        while (l < r) {
            int t = a[l];
            a[l] = a[r];
            a[r] = t;
            l++;
            r--;
        }
    }

    static void reverse(long[] a, int l, int r) {
        while (l < r) {
            long t = a[l];
            a[l] = a[r];
            a[r] = t;
            l++;
            r--;
        }
    }

    static long modPow(long b, long e, long mod) {
        long r = 1;
        b = b % mod;
        while (e > 0) {
            if ((e & 1) == 1) {
                r = (r * b) % mod;
            }
            b = (b * b) % mod;
            e >>= 1;
        }
        return r;
    }

    static class FastReader {
        BufferedReader br;
        StringTokenizer st;

        public FastReader() {
            br = new BufferedReader(new InputStreamReader(System.in));
        }

        String next() {
            while (st == null || !st.hasMoreElements()) {
                try {
                    st = new StringTokenizer(br.readLine());
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            return st.nextToken();
        }

        int nextInt() {
            return Integer.parseInt(next());
        }

        long nextLong() {
            return Long.parseLong(next());
        }

        double nextDouble() {
            return Double.parseDouble(next());
        }
    }

    static void sort_(int[] a) {
        ArrayList<Integer> ls = new ArrayList<>();
        for (int x : a)
            ls.add(x);
        Collections.sort(ls);
        for (int i = 0; i < a.length; i++)
            a[i] = ls.get(i);
    }

    static void sortRev(int[] a) {
        ArrayList<Integer> ls = new ArrayList<>();
        for (int x : a)
            ls.add(x);
        Collections.sort(ls, Collections.reverseOrder());
        for (int i = 0; i < a.length; i++)
            a[i] = ls.get(i);
    }

    static void sortRev(long[] a) {
        ArrayList<Long> ls = new ArrayList<>();
        for (long x : a)
            ls.add(x);
        Collections.sort(ls, Collections.reverseOrder());
        for (int i = 0; i < a.length; i++)
            a[i] = ls.get(i);
    }

    static <T extends Comparable<T>> void sort(ArrayList<T> list) {
        Collections.sort(list);
    }

    static <T extends Comparable<T>> void sortRev(ArrayList<T> list) {
        Collections.sort(list, Collections.reverseOrder());
    }

    static void sort_(long[] a) {
        ArrayList<Long> ls = new ArrayList<>();
        for (long x : a)
            ls.add(x);
        Collections.sort(ls);
        for (int i = 0; i < a.length; i++)
            a[i] = ls.get(i);
    }

    static void print(int[][] arr) {
        for (int[] a : arr) {
            for (int i : a)
                out.print(i + " ");
            out.println();
        }
        out.println();
    }

    static void print(long[][] arr) {
        for (long[] a : arr) {
            for (long i : a)
                out.print(i + " ");
            out.println();
        }
        out.println();
    }

    static void print(int[] a) {
        for (int i : a)
            out.print(i + " ");
        out.println();
    }

    static void print(char[] a) {
        for (char i : a)
            out.print(i + " ");
        out.println();
    }

    static void print(long[] a) {
        for (long i : a)
            out.print(i + " ");
        out.println();
    }

    static <T extends Number> void print(ArrayList<T> ls) {
        for (T i : ls)
            out.print(i + " ");
        out.println();
    }

    static int[] inputIntArr(int n) {
        int[] a = new int[n];
        for (int i = 0; i < n; i++)
            a[i] = in.nextInt();
        return a;
    }

    static long[] inputLongArr(int n) {
        long[] a = new long[n];
        for (int i = 0; i < n; i++)
            a[i] = in.nextLong();
        return a;
    }

    static int[][] inputIntArr(int n, int m) {
        int[][] a = new int[n][m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                a[i][j] = in.nextInt();
        return a;
    }

    static long[][] inputLongArr(int n, int m) {
        long[][] a = new long[n][m];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                a[i][j] = in.nextLong();
        return a;
    }

    static ArrayList<Integer> inputIntList(int n) {
        ArrayList<Integer> ls = new ArrayList<>();
        for (int i = 0; i < n; i++)
            ls.add(in.nextInt());
        return ls;
    }

    static ArrayList<Long> inputLongList(int n) {
        ArrayList<Long> ls = new ArrayList<>();
        for (int i = 0; i < n; i++)
            ls.add(in.nextLong());
        return ls;
    }
}
