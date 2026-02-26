"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import { signOut } from "next-auth/react";
import { cn } from "@/lib/utils";
import { LogOut } from "lucide-react";

const primaryNav = [
  { href: "/", label: "Сводка", icon: "\uD83C\uDFE0" },
  { href: "/messages", label: "Почта", icon: "\uD83D\uDCEC" },
];

const mainNav = [
  { href: "/directions", label: "Направления", icon: "\uD83E\uDDED" },
  { href: "/projects", label: "Проекты", icon: "\uD83D\uDCC1" },
  { href: "/experiments", label: "Эксперименты", icon: "\uD83E\uDDEA" },
  { href: "/strains", label: "Штаммы", icon: "\uD83E\uDDEC" },
  { href: "/constructs", label: "Конструкции", icon: "\uD83D\uDD27" },
  { href: "/primers", label: "Праймеры", icon: "\uD83D\uDD2C" },
  { href: "/files", label: "Файлы", icon: "\uD83D\uDCCE" },
  { href: "/protocols", label: "Протоколы", icon: "\uD83D\uDCCB" },
  { href: "/reagents", label: "Реагенты", icon: "\uD83E\uDDF4" },
  { href: "/consumables", label: "Расходники", icon: "\uD83D\uDCE6" },
  { href: "/equipment", label: "Оборудование", icon: "\u2699\uFE0F" },
  { href: "/enzymes", label: "Ферменты", icon: "\uD83D\uDD2C" },
  { href: "/team", label: "Команда", icon: "\uD83D\uDC65" },
];

interface SidebarProps {
  userName: string;
  onClose?: () => void;
}

export default function Sidebar({ userName, onClose }: SidebarProps) {
  const pathname = usePathname();

  const isActive = (href: string) => {
    if (href === "/") return pathname === "/";
    return pathname.startsWith(href);
  };

  const NavLink = ({ href, label, icon }: { href: string; label: string; icon: string }) => (
    <Link
      href={href}
      onClick={onClose}
      className={cn(
        "flex items-center gap-3 rounded-md px-3 py-2 text-sm transition-colors",
        isActive(href)
          ? "bg-slate-700 text-white"
          : "text-slate-300 hover:bg-slate-800 hover:text-white"
      )}
    >
      <span className="text-base">{icon}</span>
      <span>{label}</span>
    </Link>
  );

  return (
    <div className="flex h-full w-64 flex-col bg-slate-900">
      <div className="p-4">
        <h1 className="text-lg font-bold text-white">LIMS Lab</h1>
      </div>

      <div className="space-y-1 px-3">
        {primaryNav.map((item) => (
          <NavLink key={item.href} {...item} />
        ))}
      </div>

      <div className="my-2 border-t border-slate-700" />

      <div className="flex-1 space-y-1 overflow-y-auto px-3">
        {mainNav.map((item) => (
          <NavLink key={item.href} {...item} />
        ))}
      </div>

      <div className="border-t border-slate-700 p-4">
        <p className="mb-2 truncate text-sm text-slate-300">{userName}</p>
        <button
          onClick={() => signOut({ callbackUrl: "/login" })}
          className="flex w-full items-center gap-2 rounded-md px-3 py-2 text-sm text-slate-400 transition-colors hover:bg-slate-800 hover:text-white"
        >
          <LogOut className="h-4 w-4" />
          <span>Выйти</span>
        </button>
      </div>
    </div>
  );
}
