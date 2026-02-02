import { defineConfig } from "vite";
import path from "path";

export default defineConfig({
  // ✅ Force Vite root to be THIS folder (frontend/)
  root: __dirname,

  server: {
    port: 5173,
    strictPort: true,
    proxy: {
      "/api": {
        target: "http://127.0.0.1:8000",
        changeOrigin: true,
      },
    },
  },

  // ✅ Build output goes to project-root/dist
  build: {
    outDir: path.resolve(__dirname, "../dist"),
    emptyOutDir: true,
  },
});
