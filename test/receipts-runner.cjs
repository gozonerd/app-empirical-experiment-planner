#!/usr/bin/env node
// receipts-runner.cjs — node receipt runner (zero deps).
// Thin wrapper: loads ../src/stats-engine.js, calls StatsEngine.runReceipts(),
// prints a table, exits 1 if ANY PINNED receipt fails, exits 0 otherwise.

"use strict";

const path = require("path");
const StatsEngine = require(path.join(__dirname, "..", "src", "stats-engine.js"));

function fmt(val) {
  if (typeof val === "number") {
    if (Number.isInteger(val)) return String(val);
    return val.toFixed(6);
  }
  return String(val);
}

function pad(str, width) {
  str = String(str);
  if (str.length >= width) return str.slice(0, width - 1) + "…";
  return str + " ".repeat(width - str.length);
}

function main() {
  let rows;
  try {
    rows = StatsEngine.runReceipts();
  } catch (err) {
    console.error("FATAL: runReceipts() threw:", err && err.stack ? err.stack : err);
    process.exit(1);
    return;
  }

  const cols = [
    { key: "id", label: "ID", width: 5 },
    { key: "class", label: "CLASS", width: 7 },
    { key: "desc", label: "DESC", width: 46 },
    { key: "expected", label: "EXPECTED", width: 22 },
    { key: "got", label: "GOT", width: 20 },
    { key: "pass", label: "PASS", width: 6 },
    { key: "method", label: "METHOD", width: 40 },
    { key: "source", label: "SOURCE", width: 40 }
  ];

  const headerLine = cols.map((c) => pad(c.label, c.width)).join(" | ");
  console.log(headerLine);
  console.log("-".repeat(headerLine.length));

  let pinnedFailCount = 0;
  let pinnedCount = 0;
  let pavCount = 0;

  for (const row of rows) {
    const line = cols
      .map((c) => {
        let v = row[c.key];
        if (c.key === "got" || c.key === "expected") v = fmt(v);
        if (c.key === "pass") v = row.pass ? "PASS" : "FAIL";
        return pad(v, c.width);
      })
      .join(" | ");
    console.log(line);

    if (row.class === "PINNED") {
      pinnedCount++;
      if (!row.pass) pinnedFailCount++;
    } else if (row.class === "PAV") {
      pavCount++;
    }
  }

  console.log("-".repeat(headerLine.length));
  console.log(
    `PINNED: ${pinnedCount - pinnedFailCount}/${pinnedCount} passed. PAV (informational, not gating): ${pavCount} computed.`
  );

  if (pinnedFailCount > 0) {
    console.log(`\nFAIL: ${pinnedFailCount} pinned receipt(s) failed.`);
    process.exit(1);
  }

  console.log("\nOK: all pinned receipts pass.");
  process.exit(0);
}

main();
