import { ThemeProvider } from "../../components/theme-provider";

export async function generateStaticParams() {
  return [{ locale: "en" }, { locale: "pt" }, { locale: "es" }];
}

export default async function LocaleLayout({
  children,
  params,
}: Readonly<{
  children: React.ReactNode;
  params: Promise<{ locale: string }>;
}>) {
  await params;

  return <ThemeProvider>{children}</ThemeProvider>;
}
