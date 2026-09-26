"use client";

import { createContext, useContext, useEffect, useState } from "react";

type Theme = "light" | "dark" | "system";
type ThemeContextValue = { theme: Theme; resolvedTheme: "light" | "dark"; setTheme: (theme: Theme) => void };

const ThemeContext = createContext<ThemeContextValue | null>(null);

export function ThemeProvider({ children }: { children: React.ReactNode }) {
  const [theme, setThemeState] = useState<Theme>("system");
  const [resolvedTheme, setResolvedTheme] = useState<"light" | "dark">("dark");

  useEffect(() => {
    const saved = (localStorage.getItem("biostar-theme") as Theme | null) ?? "system";
    setThemeState(saved);
  }, []);

  useEffect(() => {
    const media = window.matchMedia("(prefers-color-scheme: dark)");
    const resolve = () => {
      const next = theme === "system" ? (media.matches ? "dark" : "light") : theme;
      setResolvedTheme(next);
      document.documentElement.dataset.theme = next;
    };
    resolve();
    media.addEventListener("change", resolve);
    return () => media.removeEventListener("change", resolve);
  }, [theme]);

  function setTheme(next: Theme) {
    setThemeState(next);
    localStorage.setItem("biostar-theme", next);
  }

  return <ThemeContext.Provider value={{ theme, resolvedTheme, setTheme }}>{children}</ThemeContext.Provider>;
}

export function useTheme() {
  const context = useContext(ThemeContext);
  if (!context) throw new Error("useTheme must be used inside ThemeProvider");
  return context;
}
