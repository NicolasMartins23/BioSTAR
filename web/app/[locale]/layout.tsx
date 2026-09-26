import { notFound } from "next/navigation";
import { ThemeProvider } from "../../components/theme-provider";

const LOCALES = ["en", "pt", "es"] as const;

export async function generateStaticParams() {
  return LOCALES.map((locale) => ({ locale }));
}

export default async function LocaleLayout({
  children,
  params,
}: Readonly<{
  children: React.ReactNode;
  params: Promise<{ locale: string }>;
}>) {
  const { locale } = await params;

  if (!LOCALES.includes(locale as (typeof LOCALES)[number])) {
    notFound();
  }

  return <ThemeProvider>{children}</ThemeProvider>;
}
