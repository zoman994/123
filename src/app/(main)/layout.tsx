import { getServerSession } from "next-auth";
import { redirect } from "next/navigation";
import { authOptions } from "@/lib/auth";
import MainLayoutClient from "./MainLayoutClient";

export default async function MainLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  const session = await getServerSession(authOptions);

  if (!session) {
    redirect("/login");
  }

  // Count unread messages for the current user
  let unreadCount = 0;
  try {
    const { prisma } = await import("@/lib/prisma");
    unreadCount = await prisma.message.count({
      where: {
        toId: session.user.id,
        status: "NEW",
      },
    });
  } catch {
    // DB might not be available yet
  }

  return (
    <MainLayoutClient
      userName={session.user.name}
      unreadCount={unreadCount}
    >
      {children}
    </MainLayoutClient>
  );
}
