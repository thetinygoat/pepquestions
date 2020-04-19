class dp {
    public static void main(String[] args) {
        solve();
    }

    public static void solve() {
        // int[][] ans = fib_logn(7);
        // int[] dp = new int[8];
        // int ans = fib_tabulation(7, dp);
        // System.out.println(ans);

        int[][] maze = { { 1, 1, 1, 1 }, { 1, 1, 0, 1 }, { 0, 1, 0, 1 }, { 1, 1, 1, 1 } };
        int[][] dp = new int[maze.length][maze[0].length];
        // int ans = mazepath_tab(maze, maze.length - 1, maze[0].length - 1, 0, 0, dp);
        // display2d(dp);
        // goldMine();
        // waysToPair();
        nInKGroups();
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
}