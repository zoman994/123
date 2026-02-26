"use client";

import { usePathname } from "next/navigation";
import Link from "next/link";
import { Menu } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";

const pageTitles: Record<string, string> = {
  "/": "Сводка",
  "/messages": "Почта",
  "/directions": "Направления",
  "/projects": "Проекты",
  "/experiments": "Эксперименты",
  "/strains": "Штаммы",
  "/constructs": "Конструкции",
  "/primers": "Праймеры",
  "/files": "Файлы",
  "/protocols": "Протоколы",
  "/reagents": "Реагенты",
  "/consumables": "Расходники",
  "/equipment": "Оборудование",
  "/enzymes": "Ферменты",
  "/team": "Команда",
};

interface TopBarProps {
  unreadCount: number;
  onMenuClick: () => void;
}

export default function TopBar({ unreadCount, onMenuClick }: TopBarProps) {
  const pathname = usePathname();

  const getTitle = () => {
    if (pageTitles[pathname]) return pageTitles[pathname];
    for (const [path, title] of Object.entries(pageTitles)) {
      if (pathname.startsWith(path) && path !== "/") return title;
    }
    return "LIMS Lab";
  };

  return (
    <header className="flex h-14 items-center justify-between border-b bg-white px-4 lg:px-6">
      <div className="flex items-center gap-3">
        <Button
          variant="ghost"
          size="icon"
          className="lg:hidden"
          onClick={onMenuClick}
        >
          <Menu className="h-5 w-5" />
        </Button>
        <h2 className="text-lg font-semibold">{getTitle()}</h2>
      </div>
      <div className="flex items-center gap-2">
        {unreadCount > 0 && (
          <Link href="/messages">
            <Badge variant="destructive" className="cursor-pointer">
              {unreadCount} новых
            </Badge>
          </Link>
        )}
      </div>
    </header>
  );
}
