import java.util.*;

class dp {
    public static void main(String[] args) {
        solve();
    }

    public static void solve() {
        // int[][] ans = fib_logn(7);
        // int[] dp = new int[8];
        // int ans = fib_tabulation(7, dp);
        // System.out.println(ans);

        // int[][] maze = { { 1, 1, 1, 1 }, { 1, 1, 0, 1 }, { 0, 1, 0, 1 }, { 1, 1, 1, 1
        // } };
        // int[][] dp = new int[15][15];
        // int ans = mazepath_tab(maze, maze.length - 1, maze[0].length - 1, 0, 0, dp);
        // display2d(dp);
        // goldMine();
        // waysToPair();
        // nInKGroups();
        // longestPalendromicSubsequenceString("BBABCBCAB");
        // longestCommonSubseqRec("geeksforgeeks", "geeksquiz", 0, 0, dp);
        // longestCommonSubseqDP();
        // editDistance();
        // stringAsSubsequence();
        // coinchnagecombination();
        int target = 100;
        int[] values = { 1, 30 };
        int[] weights = { 1, 50 };
        // int[][] dp = new int[weights.length + 1][target + 1];
        int[] dp = new int[target + 1];
        // coinChangeDP();
        // knapsack01Rec(weights, values, target, dp, 0);
        // knapsack01DP(weights, values, target, dp);
        // System.out.println(coinchangeRec(target, coins, dp, 0));
        System.out.println(unboundedKnapsackRec(weights, values, target, 0, dp));
        display1d(dp);
        // display2d(dp);
    }

    public static void display1d(int[] a) {
        for (int e : a) {
            System.out.print(e + " ");
        }
    }

    public static void display2d(int[][] arr) {
        for (int[] a : arr) {
            for (int e : a) {
                System.out.print(e + " ");
            }
            System.out.println();
        }
    }

    public static int fib_memo(int n, int[] dp) {
        if (n <= 1) {
            dp[n] = n;
            return n;
        }

        if (dp[n] != 0)
            return dp[n];

        int sum = fib_memo(n - 1, dp) + fib_memo(n - 2, dp);
        return dp[n] = sum;
    }

    public static int fib_tabulation(int n, int[] dp) {
        for (int i = 0; i <= n; i++) {
            if (i <= 1) {
                dp[i] = i;
                continue;
            }

            int sum = dp[i - 1] + dp[i - 2];
            dp[i] = sum;
        }

        return dp[n];
    }

    public static int[][] multiply_matrix(int[][] mat1, int[][] mat2) {
        int[][] res = new int[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                res[i][j] = 0;
                for (int k = 0; k < 2; k++) {
                    res[i][j] += mat1[i][k] * mat2[k][j];
                }
            }
        }

        return res;
    }

    public static int[][] fib_logn(int n) {
        if (n <= 0) {
            int[][] mat = { { 1, 1 }, { 1, 0 } };
            return mat;
        }
        int[][] mat = { { 1, 1 }, { 1, 0 } };
        if (n % 2 == 0) {
            return multiply_matrix(mat, fib_logn(n / 2));
        } else {
            return multiply_matrix(mat, multiply_matrix(mat, fib_logn(n / 2)));
        }

    }

    public static int[][] dir = { { 0, 1 }, { 1, 0 } };

    public static boolean isSafe(int[][] grid, int r, int c) {
        if (r < 0 || c < 0 || r >= grid.length || c >= grid[0].length || grid[r][c] == 0)
            return false;
        return true;
    }

    public static int mazepath(int[][] grid, int gr, int gc, int r, int c, int[][] dp) {
        if (r == gr && c == gc) {
            dp[r][c] = 1;
            return 1;
        }
        if (dp[r][c] != 0)
            return dp[r][c];
        int count = 0;
        for (int d = 0; d < dir.length; d++) {
            int nr = r + dir[d][0];
            int nc = c + dir[d][1];
            if (isSafe(grid, nr, nc)) {

                count += mazepath(grid, gr, gc, nr, nc, dp);
            }
        }
        return dp[r][c] = count;
    }

    public static int mazepath_tab(int[][] grid, int gr, int gc, int r, int c, int[][] dp) {
        for (int sr = gr; sr >= 0; sr--) {
            for (int sc = gc; sc >= 0; sc--) {
                if (sr == gr && sc == gc) {
                    dp[sr][sc] = 1;
                    continue;
                }
                int count = 0;
                for (int d = 0; d < dir.length; d++) {
                    int nr = sr + dir[d][0];
                    int nc = sc + dir[d][1];
                    if (isSafe(grid, nr, nc)) {
                        System.out.println(nr + " " + nc);
                        count += dp[nr][nc];
                    }
                }
                dp[sr][sc] = count;
            }
        }
        return dp[0][0];
    }

    int[] moves = { 1, 2, 3, 4, 5, 6 };
    public static int[][] mine = { { 1, 3, 1, 5 }, { 2, 2, 4, 1 }, { 5, 0, 2, 3 }, { 0, 6, 1, 2 } };

    public static int goldMine_(int r, int c, int m, int n, int[][] dp) {
        if (c == n - 1) {
            dp[r][c] = mine[r][c];
            return dp[r][c];
        }
        if (dp[r][c] != 0)
            return dp[r][c];
        int right = 0, bottom = 0, diag = 0;
        if (c + 1 < n)
            right = goldMine_(r, c + 1, m, n, dp);
        if (r + 1 < m && c + 1 < n)
            diag = goldMine_(r + 1, c + 1, m, n, dp);
        if (r - 1 >= 0 && c + 1 < n)
            bottom = goldMine_(r - 1, c + 1, m, n, dp);

        return dp[r][c] = Math.max(right, Math.max(bottom, diag)) + mine[r][c];
    }

    public static int goldMineDp(int m, int n, int[][] dp) {
        for (int c = n - 1; c >= 0; c--) {
            for (int r = m - 1; r >= 0; r--) {
                if (c == n - 1 && r == m - 1) {
                    dp[r][c] = mine[r][c];
                    continue;
                }
                int right = 0, bottom = 0, diag = 0;
                if (c + 1 < n)
                    right = dp[r][c + 1];
                if (r + 1 < m && c + 1 < n)
                    diag = dp[r + 1][c + 1];
                if (r - 1 >= 0 && c + 1 < n)
                    bottom = dp[r - 1][c + 1];

                dp[r][c] = Math.max(right, Math.max(bottom, diag)) + mine[r][c];
            }
        }
        return dp[0][0];
    }

    public static void goldMine() {
        int m = mine.length;
        int n = mine[0].length;
        int[][] dp = new int[m][n];
        int max = 0;
        // for (int i = 0; i < m; i++) {
        // max = Math.max(max, goldMine_(i, 0, m, n, dp));
        // }
        goldMineDp(m, n, dp);
        display2d(dp);
        System.out.println(max);
    }

    public static int waysToPair_(int c, int[] dp) {
        if (c <= 1) {
            return dp[c] = 1;
        }
        if (dp[c] != 0)
            return dp[c];
        int single = waysToPair_(c - 1, dp);
        int pair = waysToPair_(c - 2, dp) * c - 1;
        return dp[c] = single + pair;
    }

    public static void waysToPair() {
        int p = 3;
        int[] dp = new int[p + 1];
        waysToPair_(p, dp);
        display1d(dp);
    }

    public static int nInKGroups_(int n, int k, int[][] dp) {
        if (n < k)
            return 0;
        if (k == 1 || n == k) {
            return dp[n][k] = 1;
        }
        if (dp[n][k] != 0)
            return dp[n][k];
        int count = 0;
        count += nInKGroups_(n - 1, k - 1, dp);
        count += nInKGroups_(n - 1, k, dp) * k;
        return dp[n][k] = count;
    }

    public static void nInKGroups() {
        int n = 3;
        int k = 2;
        int[][] dp = new int[n + 1][k + 1];
        int res = nInKGroups_(n, k, dp);
        display2d(dp);
        System.out.println(res);
    }

    // ================================= Strings
    // =====================================

    public static int longestPalendromicSubstring(String s) {
        int[][] dp = new int[s.length()][s.length()];
        int maxlen = 0;
        for (int gap = 0; gap < s.length(); gap++) {
            for (int si = 0, ei = gap; ei < s.length(); si++, ei++) {
                if (gap == 0)
                    dp[si][ei] = 1;
                else if (gap == 1 && s.charAt(si) == s.charAt(ei))
                    dp[si][ei] = 2;
                else if (s.charAt(si) == s.charAt(ei) && dp[si + 1][ei - 1] != 0) {
                    dp[si][ei] = dp[si + 1][ei - 1] + 2;
                    maxlen = Math.max(maxlen, dp[si][ei]);
                }
            }
        }
        return maxlen;
    }

    // public static int longestPalendromicSubsequence(String s) {
    // int dp[][] = new int[s.length()][s.length()];
    // }
    public static void longestPalendromicSubsequenceString(String str) {
        StringBuilder dp[][] = new StringBuilder[str.length()][str.length()];
        for (int i = 0; i < dp.length; i++) {
            for (int j = 0; j < dp[0].length; j++) {
                dp[i][j] = new StringBuilder("");
            }
        }
        for (int gap = 0; gap < str.length(); gap++) {
            for (int si = 0, ei = gap; ei < str.length(); ei++, si++) {
                if (gap == 0)
                    dp[si][ei].append(str.charAt(si));
                else if (gap == 1 && str.charAt(si) == str.charAt(ei))
                    dp[si][ei].append(str.substring(si, ei + 1));
                else if (str.charAt(si) == str.charAt(ei)) {
                    dp[si][ei].append(str.charAt(si));
                    dp[si][ei].append(dp[si + 1][ei - 1]);
                    dp[si][ei].append(str.charAt(ei));
                } else {
                    dp[si][ei].append(
                            dp[si + 1][ei].length() >= dp[si][ei - 1].length() ? dp[si + 1][ei] : dp[si][ei - 1]);
                }
            }
        }

        System.out.println(dp[0][str.length() - 1].toString());
    }

    public static int longestCommonSubseqRec(String s1, String s2, int i, int j, int[][] dp) {
        if (i == s1.length() || j == s2.length())
            return 0;
        if (dp[i][j] != 0)
            return dp[i][j];
        if (s1.charAt(i) == s2.charAt(j))
            return dp[i][j] = longestCommonSubseqRec(s1, s2, i + 1, j + 1, dp) + 1;
        int l = longestCommonSubseqRec(s1, s2, i + 1, j, dp);
        int r = longestCommonSubseqRec(s1, s2, i, j + 1, dp);

        return dp[i][j] = Math.max(l, r);
    }

    public static int longestCommonSubseqDP() {
        String s1 = "geeksforgeeks";
        String s2 = "geeksquiz";
        int[][] dp = new int[s1.length()][s2.length()];
        for (int i = s1.length() - 1; i >= 0; i--) {
            for (int j = s2.length() - 1; j >= 0; j--) {
                if (s1.charAt(i) == s2.charAt(j))
                    dp[i][j] = dp[i + 1][j + 1] + 1;
                else if (s1.charAt(i) != s2.charAt(j)) {
                    int l = dp[i + 1][j];
                    int r = dp[i][j + 1];
                    dp[i][j] = Math.max(l, r);
                }
            }
        }
        display2d(dp);
        return dp[0][0];
    }

    public static int max = 0;

    public static int longestCommonSubstringRec(String s1, String s2, int i, int j, int[][] dp) {
        if (i == s1.length() || j == s2.length())
            return 0;
        int a = 0;
        if (s1.charAt(i) == s2.charAt(j)) {
            a = longestCommonSubstringRec(s1, s2, i + 1, j + 1, dp) + 1;
            max = Math.max(a, max);
        }
        longestCommonSubstringRec(s1, s2, i + 1, j, dp);
        longestCommonSubstringRec(s1, s2, i, j + 1, dp);
        return dp[i][j] = a;
    }

    public static int longestCommonSubstringDP(String s1, String s2, int i, int j, int[][] dp) {
        int max = 0;
        for (i = s1.length() - 1; i >= 0; i--) {
            for (j = s2.length() - 1; j >= 0; j--) {
                if (s1.charAt(i) == s2.charAt(j)) {
                    dp[i][j] = dp[i + 1][j + 1] + 1;
                    max = Math.max(dp[i][j], max);
                }
            }
        }
        return dp[0][0];
    }

    public static void editDistance() {
        String s1 = "horse";
        String s2 = "ros";

        int[][] dp = new int[s1.length() + 1][s2.length() + 1];

        for (int i = 0; i < dp.length; i++) {
            for (int j = 0; j < dp[0].length; j++) {
                if (i == 0)
                    dp[i][j] = j;
                else if (j == 0)
                    dp[i][j] = i;
                else if (s1.charAt(i - 1) == s2.charAt(j - 1))
                    dp[i][j] = dp[i - 1][j - 1];
                else
                    dp[i][j] = Math.min(dp[i - 1][j - 1], Math.min(dp[i][j - 1], dp[i - 1][j])) + 1;
            }
        }

        display2d(dp);

    }

    public static void stringAsSubsequence() {
        String s1 = "geeksforgeeks";
        String s2 = "gks";
        int[][] dp = new int[s1.length()][s2.length()];
        for (int i = s1.length() - 1; i >= 0; i--) {
            for (int j = s2.length() - 1; j >= 0; j--) {
                if (i == s1.length() - 1 && j == s2.length() - 1)
                    dp[i][j] = 1;
                else if (i == s1.length() - 1 && j != s2.length() - 1)
                    dp[i][j] = 0;
                else if (i != s1.length() - 1 && j == s2.length() - 1)
                    dp[i][j] = 0;
                else if (s1.charAt(i) == s2.charAt(j))
                    dp[i][j] = dp[i + 1][j + 1] + dp[i + 1][j];
            }
        }
        display2d(dp);
    }

    // ============= coin change ================

    public static void coinchnageperm() {
        int target = 10;
        int[] coins = { 2, 3, 5, 7 };
        int[] dp = new int[target + 1];
        for (int i = 0; i < dp.length; i++) {
            if (i == 0) {
                dp[i] = 1;
                continue;
            }
            for (int c = 0; c < coins.length; c++) {
                if (i - coins[c] >= 0) {
                    dp[i] += dp[i - coins[c]];
                }
            }
        }

        display1d(dp);
    }

    public static void coinchnagecombination() {
        int target = 10;
        int[] coins = { 2, 3, 5, 7 };
        int[] dp = new int[target + 1];
        dp[0] = 1;
        for (int c = 0; c < coins.length; c++) {
            for (int i = 0; i < dp.length; i++) {
                if (i - coins[c] >= 0) {
                    dp[i] += dp[i - coins[c]];
                }
            }
        }

        display1d(dp);
    }

    // =============== solution of linear equation
    public static void lineareq() {
        int tar = 10;
        int[] coeff = { 2, 3, 5, 7 };
        int[] dp = new int[tar + 1];
        dp[0] = 1;
        for (int c : coeff) {
            for (int t = 0; t <= tar; t++) {
                if (t - c >= 0) {
                    dp[t] += dp[t - c];
                }
            }
        }
    }

    // ==================== coin change with single usage

    public static int coinchangeRec(int t, int[] coins, int[][] dp, int idx) {
        if (t == 0 || idx == coins.length) {
            return dp[idx][t] = t == 0 ? 1 : 0;
        }

        if (dp[idx][t] != 0)
            return dp[idx][t];
        int count = 0;
        if (t - coins[idx] >= 0) {
            count += coinchangeRec(t - coins[idx], coins, dp, idx + 1);
        }
        count += coinchangeRec(t, coins, dp, idx + 1);
        return dp[idx][t] = count;
    }

    public static void coinChangeDP() {
        int t = 10;
        int[] coins = { 2, 3, 1, 5, 6 };
        int[][] dp = new int[coins.length + 1][t + 1];
        dp[0][0] = 1;

        for (int i = 1; i < dp.length; i++) {
            for (int j = 0; j <= t; j++) {
                int count = 0;
                if (j - coins[i - 1] >= 0) { // i-1 because we add 0 row
                    count += dp[i - 1][j - coins[i - 1]];
                }

                count += dp[i - 1][j];
                dp[i][j] = count;
            }
        }

        display2d(dp);
    }

    // =============== knapsack 0-1

    public static int knapsack01Rec(int[] weights, int[] values, int W, int[][] dp, int idx) {
        if (W == 0 || idx == weights.length) {
            return dp[idx][W] = 0;
        }

        if (dp[idx][W] != 0)
            return dp[idx][W];

        int t = 0, nt = 0;
        if (W - weights[idx] >= 0) {
            t = knapsack01Rec(weights, values, W - weights[idx], dp, idx + 1) + values[idx];
        }
        nt = knapsack01Rec(weights, values, W, dp, idx + 1);

        return dp[idx][W] = Math.max(t, nt);
    }

    public static void knapsack01DP(int[] weights, int[] values, int W, int[][] dp) {
        for (int i = dp.length - 1; i >= 0; i--) {
            for (int j = 0; j <= W; j++) {
                if (j == 0 || i == dp.length - 1) {
                    dp[i][j] = 0;
                    continue;
                }

                int t = 0, nt = 0;
                if (j - weights[i] >= 0) {
                    t = dp[i + 1][j - weights[i]] + values[i];
                }
                nt = dp[i + 1][j];
                dp[i][j] = Math.max(t, nt);
            }
        }
    }
    // ================ unbounded knapsack

    public static int unboundedKnapsackRec(int[] weights, int[] values, int W, int idx, int[] dp) {
        if (W == 0 || idx == weights.length) {
            return dp[idx] = 1;
        }

        if (dp[idx] != 0)
            return dp[idx];
        int count = 0;
        for (int w : weights) {
            if (W - w >= 0) {
                count += unboundedKnapsackRec(weights, values, W - w, idx + 1, dp) + values[idx];
            }
        }

        return dp[idx] = count;
    }

}