import "../globals.css";

export async function generateStaticParams() {
  return [{ locale: "en" }, { locale: "pt" }, { locale: "es" }];
}

export default async function LocaleLayout({ children }: Readonly<{ children: React.ReactNode; params: Promise<{ locale: string }> }>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
