# GenoScope Frontend

## Установка и запуск

### 1. Установка зависимостей
```bash
cd /home/asd/.virtualenvs/BioForge_edit_branch/frontend
npm install
```

### 2. Запуск приложения
```bash
npm start
```

Приложение будет доступно по адресу: http://localhost:3000

### 3. Если возникнут ошибки с зависимостями:
```bash
# Очистка кэша
rm -rf node_modules package-lock.json
npm cache clean --force

# Переустановка
npm install
```

### 4. Тестовый вход
- Email: demo@genoscope.com
- Password: demo123

Или используйте кнопку "Use Demo Account" на странице входа.

## Структура проекта

```
frontend/
├── public/
│   ├── index.html
│   └── manifest.json
├── src/
│   ├── components/
│   │   └── Layout/
│   │       ├── Header.tsx
│   │       └── Sidebar.tsx
│   ├── pages/
│   │   ├── Login.tsx
│   │   ├── Dashboard.tsx
│   │   ├── Analysis.tsx
│   │   └── ... (другие страницы)
│   ├── App.tsx
│   ├── index.tsx
│   └── index.css
├── package.json
└── tsconfig.json
```

## Возможные проблемы и решения

### Ошибка "Module not found"
Убедитесь, что все зависимости установлены:
```bash
npm install
```

### Ошибка с портом 3000
Если порт занят:
```bash
PORT=3001 npm start
```

### Ошибки TypeScript
Проверьте версию TypeScript:
```bash
npm list typescript
```
Должна быть версия 4.9.5 или выше.
