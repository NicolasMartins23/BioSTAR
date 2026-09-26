"use client";

import { useTheme } from "./theme-provider";

export function ThemeToggle() {
  const { theme, setTheme } = useTheme();
  return (
    <div className="theme-toggle" aria-label="Theme">
      {([['light','☼'], ['system','◐'], ['dark','◑']] as const).map(([value, icon]) => (
        <button key={value} className={theme === value ? "active" : ""} onClick={() => setTheme(value)} title={value} aria-label={value}>
          {icon}
        </button>
      ))}
    </div>
  );
}
