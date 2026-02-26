import type {
  Organism,
  DirectionType,
  ProjectStatus,
  ExpStatus,
  ExpResult,
  WorkType,
  Priority,
  StrainType,
  ConstructType,
  GelType,
  AnalysisType,
  StockType,
  SeqStatus,
  Verification,
  ProtocolCategory,
  FileType,
  MessageStatus,
  MessageColor,
  SampleType,
  ReagentCategory,
  ConsumableCategory,
} from "@/generated/prisma/enums";

export const ORGANISM_LABELS: Record<Organism, string> = {
  AN: "A. niger",
  TR: "T. reesei",
  KL: "K. lactis",
  SC: "S. cerevisiae",
  PP: "P. pastoris",
  YL: "Y. lipolytica",
};

export const DIRECTION_TYPE_LABELS: Record<DirectionType, string> = {
  FUNDAMENTAL: "Фундаментальное",
  APPLIED: "Прикладное",
};

export const PROJECT_STATUS_LABELS: Record<ProjectStatus, string> = {
  PLANNED: "Планируется",
  IN_PROGRESS: "В работе",
  PAUSED: "Приостановлен",
  COMPLETED: "Завершён",
  CANCELLED: "Отменён",
};

export const EXP_STATUS_LABELS: Record<ExpStatus, string> = {
  PLANNED: "Запланирован",
  IN_PROGRESS: "В работе",
  COMPLETED: "Завершён",
  CANCELLED: "Отменён",
};

export const EXP_RESULT_LABELS: Record<ExpResult, string> = {
  SUCCESS: "Успешно",
  PARTIAL: "Частично",
  FAILURE: "Неудачно",
  UNDEFINED: "Не определён",
};

export const WORK_TYPE_LABELS: Record<WorkType, string> = {
  CLONING: "Клонирование",
  TRANSFORMATION: "Трансформация",
  FERMENTATION: "Ферментация",
  ANALYTICS: "Аналитика",
  SCREENING: "Скрининг",
};

export const PRIORITY_LABELS: Record<Priority, string> = {
  HIGH: "Высокий",
  MEDIUM: "Средний",
  LOW: "Низкий",
};

export const STRAIN_TYPE_LABELS: Record<StrainType, string> = {
  WILD_TYPE: "Дикий тип",
  MUTANT: "Мутант",
  RECOMBINANT: "Рекомбинант",
  MUSEUM: "Музейный",
};

export const CONSTRUCT_TYPE_LABELS: Record<ConstructType, string> = {
  EXPRESSION_VECTOR: "Экспрессионный вектор",
  KNOCKOUT_CASSETTE: "Нокаут-кассета",
  INTEGRATION_PLASMID: "Интеграционная плазмида",
};

export const GEL_TYPE_LABELS: Record<GelType, string> = {
  AGAROSE: "Агарозный",
  PAAG: "ПААГ",
  SDS_PAGE: "SDS-PAGE",
  NATIVE: "Нативный",
};

export const ANALYSIS_TYPE_LABELS: Record<AnalysisType, string> = {
  DNA: "ДНК",
  RNA: "РНК",
  PROT: "Белок",
};

export const STOCK_TYPE_LABELS: Record<StockType, string> = {
  GLYCEROL: "Глицериновый сток",
  SPORES: "Споры",
  MYCELIUM: "Мицелий",
  CRYOVIAL: "Криопробирка",
};

export const SEQ_STATUS_LABELS: Record<SeqStatus, string> = {
  NOT_SEQUENCED: "Не секвенирован",
  IN_PROGRESS: "В работе",
  CONFIRMED: "Подтверждён",
  ERROR: "Ошибка",
};

export const VERIFICATION_LABELS: Record<Verification, string> = {
  PCR: "ПЦР",
  SEQUENCING: "Секвенирование",
  RESTRICTION: "Рестрикция",
  NOT_VERIFIED: "Не проверен",
};

export const PROTOCOL_CATEGORY_LABELS: Record<ProtocolCategory, string> = {
  CLONING: "Клонирование",
  TRANSFORMATION: "Трансформация",
  FERMENTATION: "Ферментация",
  ANALYTICS: "Аналитика",
};

export const FILE_TYPE_LABELS: Record<FileType, string> = {
  GEL: "Гель",
  SEQUENCE: "Последовательность",
  DATA: "Данные",
  REPORT: "Отчёт",
  PHOTO: "Фото",
  MAP: "Карта",
  PROTOCOL: "Протокол",
  OTHER: "Другое",
};

export const MESSAGE_STATUS_LABELS: Record<MessageStatus, string> = {
  NEW: "Новое",
  READ: "Прочитано",
  IN_PROGRESS: "В работе",
  DONE: "Выполнено",
};

export const MESSAGE_COLOR_LABELS: Record<MessageColor, string> = {
  WHITE: "Белый",
  YELLOW: "Жёлтый",
  RED: "Красный",
  GREEN: "Зелёный",
};

export const SAMPLE_TYPE_LABELS: Record<SampleType, string> = {
  GENOMIC_DNA: "Геномная ДНК",
  PLASMID: "Плазмида",
  PCR_PRODUCT: "ПЦР-продукт",
  CDNA: "кДНК",
  RNA: "РНК",
};

export const REAGENT_CATEGORY_LABELS: Record<ReagentCategory, string> = {
  ENZYMES: "Ферменты",
  BUFFERS: "Буферы",
  ANTIBIOTICS: "Антибиотики",
  DYES: "Красители",
  MEDIA: "Среды",
};

export const CONSUMABLE_CATEGORY_LABELS: Record<ConsumableCategory, string> = {
  PLASTIC: "Пластик",
  FILTERS: "Фильтры",
  GLOVES: "Перчатки",
  TIPS: "Наконечники",
};

// Цвета для статусных бейджей
export const STATUS_COLORS: Record<string, string> = {
  // ProjectStatus / ExpStatus
  "Планируется": "gray",
  "Запланирован": "gray",
  "В работе": "yellow",
  "Приостановлен": "orange",
  "Завершён": "green",
  "Отменён": "red",
  // ExpResult
  "Успешно": "green",
  "Частично": "blue",
  "Неудачно": "red",
  "Не определён": "gray",
  // SeqStatus
  "Не секвенирован": "gray",
  "Подтверждён": "green",
  "Ошибка": "red",
  // MessageStatus
  "Новое": "blue",
  "Прочитано": "gray",
  "Выполнено": "green",
  // Equipment / general
  "Работает": "green",
  "На ремонте": "yellow",
  "Списан": "red",
  "Активный": "green",
  "Удалён": "red",
  // Reagent
  "В наличии": "green",
  "Мало": "yellow",
  "Закончился": "red",
};
