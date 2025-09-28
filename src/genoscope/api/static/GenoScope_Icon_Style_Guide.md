# GenoScope — Icon Style Guide (v1)

## Стиль
- **Outline** 1.75px, `round` cap/join, монохром `currentColor` для UI.
- **Solid-версии** для микро-размеров (16–20 px) — по мере необходимости.
- **Duotone** — только для маркетинга/пустых состояний.

## Сетка и геометрия
- Базовая сетка: **24×24** px.
- Общие параметры: единый радиус (оптически ~2 px), диагонали под 45°.
- Допускаются оптические «вылеты» за сетку до 0.5 px.

## Палитра (рекомендация)
- Базовый цвет иконок: `#0F172A` (или `currentColor`).
- Акцент для промо: градиент `#22D3EE → #6366F1`.

## Экспорт
- SVG (основа), дополнительно PNG 16/20/24/32/48 при необходимости.
- Favicon: `favicon.svg` (есть в комплекте). PWA: `site.webmanifest` (шаблон в комплекте).

## Состояния
- Default: 70% непрозрачность, Hover: 100%, Disabled: 40%.
- Темная тема: следить за контрастом ≥ 4.5:1 для ключевых иконок.

## Набор v1 (16 штук)
- Доменные: dna-scope, dna, reads, alignment, crispr, protein, chromosome, pipette, flask.
- Системные: home, search, settings, upload, download, filter, table.

## Встраивание
```html
<link rel="icon" href="/favicon.svg" type="image/svg+xml">
<link rel="manifest" href="/site.webmanifest">
<meta name="theme-color" content="#0F172A">

<!-- Пример иконки как React-компонент -->
<svg class="w-6 h-6 text-slate-800">{/* вставь содержимое нужного SVG */}</svg>
```

## Замечания по доступности
- Не полагаться только на цвет для смысла.
- Предусмотреть подписи/aria-label для action-иконок.

---

© GenoScope Icons v1