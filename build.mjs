// Assembly: inline src/stats-engine.js into src/app-shell.html at the marker → index.html
import { readFileSync, writeFileSync } from "node:fs";

const MARKER = "/*__INLINE_STATS_ENGINE__*/";
const shell = readFileSync("src/app-shell.html", "utf8");
const engine = readFileSync("src/stats-engine.js", "utf8");

if (!shell.includes(MARKER)) {
  console.error("build failed: marker not found in src/app-shell.html");
  process.exit(1);
}
writeFileSync("index.html", shell.replace(MARKER, () => engine));
console.log("built index.html (" + Math.round((shell.length + engine.length) / 1024) + " KB)");
