import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.abs;
import java.io.*;
import java.util.*;

public class Main {
    public static PrintWriter out = new PrintWriter(new BufferedOutputStream(System.out));
    static FastReader in = new FastReader();
    // static final int MOD = (int) 1e9 + 7;
    static final int MOD = 998244353;
    static List<List<int[]>> adj;
    private static final int MAX_LOG = 13;
    private static final int N = (int) 1e5 + 1;
    static int[][] par = new int[N][MAX_LOG];
    static int[] dist = new int[N];
    static int[] depth = new int[N];

    public static void main(String[] Hi) {
        int T = in.nextInt();
        while (T-- > 0) {
            solve();
        }
        out.close();
    }

    static void solve() {
        int n = in.nextInt();
        adj = new ArrayList<>();
        for (int i = 0; i <= n; i++) {
            adj.add(new ArrayList<>());
        }
        for (int i = 0; i < n - 1; i++) {
            int a = in.nextInt();
            int b = in.nextInt();
            int c = in.nextInt();
            adj.get(a).add(new int[] { b, c });
            adj.get(b).add(new int[] { a, c });
        }
        dfs(1, 0);
        while (true) {
            String str = in.next();
            if (str.equals("DONE")) {
                break;
            }
            if (str.equals("DIST")) {
                int a = in.nextInt();
                int b = in.nextInt();
                int ans = dist[a] + dist[b];
                int lca = LCA_(a, b);
                ans -= 2 * (dist[lca]);
                out.println(ans);
            } else {
                int a = in.nextInt();
                int b = in.nextInt();
                int k = in.nextInt();
                if (k == 1) {
                    out.println(a);
                    continue;
                }
                int lca = LCA_(a, b);
                int d = depth[a] + depth[b] - (2 * depth[LCA_(a, b)]);
                int distA_lca = depth[a] - depth[lca];
                int distB_lca = depth[b] - depth[lca];
                k--;
                int ans = 0;
                if (k <= distA_lca) {
                    ans = Kthparent(a, k);
                } else {
                    int newK = distB_lca - (k - distA_lca);
                    ans = Kthparent(b, newK);
                }

                out.println(ans);
            }

        }
    }

    static int Kthparent(int node, int k) {
        for (int i = MAX_LOG - 1; i >= 0; i--) {
            if ((k & (1 << i)) != 0) {
                node = par[node][i];
                if (node == 0)
                    break;
            }
        }
        return node;
    }

    static void dfs(int node, int parent) {
        par[node][0] = parent;
        depth[node] = 1 + depth[parent];
        for (int j = 1; j < MAX_LOG; j++) {
            par[node][j] = par[par[node][j - 1]][j - 1];
        }
        for (int[] a : adj.get(node)) {
            int adjNode = a[0];
            int distance = a[1];
            if (adjNode != parent) {
                dist[adjNode] = dist[node] + distance;
                dfs(adjNode, node);
            }
        }
    }

    private static int LCA_(int nodeA, int nodeB) {
        if (nodeA == nodeB)
            return nodeA;

        if (depth[nodeA] < depth[nodeB]) {
            int tempNode = nodeA;
            nodeA = nodeB;
            nodeB = tempNode;
        }

        int nodeDiff = depth[nodeA] - depth[nodeB];
        for (int j = MAX_LOG - 1; j >= 0; j--) {
            if (((1 << j) & nodeDiff) != 0) {
                nodeA = par[nodeA][j];
            }
        }

        for (int j = MAX_LOG - 1; j >= 0; j--) {
            if (par[nodeA][j] != par[nodeB][j]) {
                nodeA = par[nodeA][j];
                nodeB = par[nodeB][j];
            }
        }

        return (nodeA != nodeB ? par[nodeA][0] : nodeA);
    }

    /*-------------------------------------------------------------------------------------------------------------- */

    class Pair {
        int first, second;

        Pair(int f, int s) {
            first = f;
            second = s;
        }

    }

    static long gcd(long a, long b) {
        if (a == 0)
            return b;
        return gcd(b % a, a);
    }

    static long lcm(long a, long b) {
        return Math.abs(a * b) / gcd(a, b);
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

    static void sort(int[] a) {
        ArrayList<Integer> ls = new ArrayList<Integer>();
        for (int x : a)
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
    }

    static void print(long[][] arr) {
        for (long[] a : arr) {
            for (long i : a)
                out.print(i + " ");
            out.println();
        }
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
