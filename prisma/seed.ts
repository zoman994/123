import { PrismaClient } from "../src/generated/prisma/client";
import { PrismaPg } from "@prisma/adapter-pg";
import { Pool } from "pg";
import { hashSync } from "bcryptjs";

import "dotenv/config";

const pool = new Pool({ connectionString: process.env.DATABASE_URL });
const adapter = new PrismaPg(pool);
const prisma = new PrismaClient({ adapter });

async function main() {
  console.log("Seeding database...");

  // Tags
  const tagNames = [
    "GLA", "BGL", "CHY", "PYR",
    "секреция", "промотор", "оптимизация", "морфология", "масштабирование",
  ];
  for (const name of tagNames) {
    await prisma.tag.upsert({
      where: { name },
      update: {},
      create: { name },
    });
  }
  console.log("Tags seeded");

  // Storage Locations
  const locationNames = ["Холодильник 1", "Морозилка -20", "Морозилка -80", "Шкаф"];
  for (const name of locationNames) {
    await prisma.storageLocation.upsert({
      where: { name },
      update: {},
      create: { name },
    });
  }
  console.log("Storage locations seeded");

  // Boxes
  for (let i = 1; i <= 10; i++) {
    const name = `Бокс-${i}`;
    await prisma.box.upsert({
      where: { name },
      update: {},
      create: { name },
    });
  }
  console.log("Boxes seeded");

  // Promoters
  const promoterNames = ["PglaA", "PgpdA", "Ptef1", "PalcA"];
  for (const name of promoterNames) {
    await prisma.promoter.upsert({
      where: { name },
      update: {},
      create: { name },
    });
  }
  console.log("Promoters seeded");

  // Terminators
  const terminatorNames = ["TtrpC", "TglaA", "Tcyc1"];
  for (const name of terminatorNames) {
    await prisma.terminator.upsert({
      where: { name },
      update: {},
      create: { name },
    });
  }
  console.log("Terminators seeded");

  // Genes
  const geneNames = ["glaA", "bgl1", "cbh1", "pyrG"];
  for (const name of geneNames) {
    await prisma.gene.upsert({
      where: { name },
      update: {},
      create: { name },
    });
  }
  console.log("Genes seeded");

  // Admin user
  const adminEmail = "admin@lab.local";
  const passwordHash = hashSync("admin123", 10);
  await prisma.user.upsert({
    where: { email: adminEmail },
    update: {},
    create: {
      email: adminEmail,
      passwordHash,
      fullName: "Администратор",
      role: "Администратор",
      status: "Активный",
    },
  });
  console.log("Admin user seeded");

  console.log("Database seeding completed!");
}

main()
  .catch((e) => {
    console.error(e);
    process.exit(1);
  })
  .finally(async () => {
    await prisma.$disconnect();
  });
