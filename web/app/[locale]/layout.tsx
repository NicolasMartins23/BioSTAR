import "../globals.css";
import { ThemeProvider } from "../../components/theme-provider";

export async function generateStaticParams() {
  return [{ locale: "en" }, { locale: "pt" }, { locale: "es" }];
}

export default async function LocaleLayout({ children, params }: Readonly<{ children: React.ReactNode; params: Promise<{ locale: string }> }>) {
  const { locale } = await params;
  const safeLocale = ["en", "pt", "es"].includes(locale) ? locale : "en";
  return <html lang={safeLocale} suppressHydrationWarning><body><ThemeProvider>{children}</ThemeProvider></body></html>;
}
