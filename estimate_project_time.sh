#!/bin/bash
# Script to estimate time spent on this specific project in a shared Git repository
# Usage: ./estimate_project_time.sh [author_name]

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get author name from argument or git config
if [ -n "$1" ]; then
    AUTHOR="$1"
else
    AUTHOR=$(git config user.name)
fi

PROJECT_DIR="."

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}   Project Time Estimation Tool${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo -e "${GREEN}Author:${NC} $AUTHOR"
echo -e "${GREEN}Directory:${NC} $(pwd)"
echo -e "${GREEN}Analysis Date:${NC} $(date +%Y-%m-%d)"
echo ""

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    echo "Error: Not in a git repository"
    exit 1
fi

# Check if author has any commits
TOTAL_COMMITS=$(git log --author="$AUTHOR" --oneline -- $PROJECT_DIR 2>/dev/null | wc -l)
if [ "$TOTAL_COMMITS" -eq 0 ]; then
    echo -e "${YELLOW}Warning: No commits found for author '$AUTHOR' in this directory${NC}"
    echo "Try running: ./estimate_project_time.sh \"Your Name\""
    exit 1
fi

echo -e "${BLUE}--- Commit Timeline ---${NC}"

# First and last commits
FIRST=$(git log --author="$AUTHOR" --reverse --pretty=format:"%ad" --date=short -- $PROJECT_DIR | head -1)
LAST=$(git log --author="$AUTHOR" --pretty=format:"%ad" --date=short -- $PROJECT_DIR | head -1)

echo -e "${GREEN}First commit:${NC} $FIRST"
echo -e "${GREEN}Last commit:${NC} $LAST"

# Calculate calendar days span
if [[ -n "$FIRST" && -n "$LAST" ]]; then
    FIRST_EPOCH=$(date -d "$FIRST" +%s 2>/dev/null || date -j -f "%Y-%m-%d" "$FIRST" +%s)
    LAST_EPOCH=$(date -d "$LAST" +%s 2>/dev/null || date -j -f "%Y-%m-%d" "$LAST" +%s)
    CALENDAR_DAYS=$(( (LAST_EPOCH - FIRST_EPOCH) / 86400 ))
    echo -e "${GREEN}Calendar span:${NC} $CALENDAR_DAYS days"
fi

echo ""
echo -e "${BLUE}--- Activity Metrics ---${NC}"

# Active working days (unique commit dates)
ACTIVE_DAYS=$(git log --author="$AUTHOR" --pretty=format:"%ad" --date=short -- $PROJECT_DIR | sort | uniq | wc -l)
echo -e "${GREEN}Active working days:${NC} $ACTIVE_DAYS days"

# Total commits
echo -e "${GREEN}Total commits:${NC} $TOTAL_COMMITS"

# Average commits per active day
if [ "$ACTIVE_DAYS" -gt 0 ]; then
    AVG_COMMITS=$(echo "scale=1; $TOTAL_COMMITS / $ACTIVE_DAYS" | bc)
    echo -e "${GREEN}Avg commits/day:${NC} $AVG_COMMITS"
fi

echo ""
echo -e "${BLUE}--- Time Estimates ---${NC}"

# Calculate estimated hours with different assumptions
HOURS_LOW=$((ACTIVE_DAYS * 2))
HOURS_MED=$((ACTIVE_DAYS * 3))
HOURS_HIGH=$((ACTIVE_DAYS * 4))

echo -e "${GREEN}Estimated hours:${NC}"
echo "  Conservative (2hrs/day): $HOURS_LOW hours"
echo "  Moderate (3hrs/day):     $HOURS_MED hours"
echo "  Intensive (4hrs/day):    $HOURS_HIGH hours"

echo ""
echo -e "${BLUE}--- Monthly Activity ---${NC}"
git log --author="$AUTHOR" --pretty=format:"%ad" --date=format:'%Y-%m' -- $PROJECT_DIR | \
    sort | uniq -c | \
    awk '{printf "  %s: %2d days\n", $2, $1}'

echo ""
echo -e "${BLUE}--- Recent Commits (Last 10) ---${NC}"
git log --author="$AUTHOR" --pretty=format:"%ad %s" --date=short -- $PROJECT_DIR | head -10

echo ""
echo -e "${BLUE}--- File Type Breakdown ---${NC}"

# Count commits by file type
echo "Commits by file extension:"
git log --author="$AUTHOR" --name-only --pretty=format: -- $PROJECT_DIR | \
    grep -v '^$' | \
    sed 's/.*\.//' | \
    sort | uniq -c | sort -rn | head -10 | \
    awk '{printf "  %s: %d commits\n", $2, $1}'

echo ""
echo -e "${BLUE}--- Summary ---${NC}"
echo -e "${GREEN}Project:${NC} Parse Single-Cell RNA-seq Analysis"
echo -e "${GREEN}Period:${NC} $FIRST to $LAST"
echo -e "${GREEN}Estimated Total Time:${NC} ${HOURS_LOW}-${HOURS_HIGH} hours (${ACTIVE_DAYS} active days)"

# Calculate work weeks (assuming 20 hours/week for research projects)
WEEKS_LOW=$(echo "scale=1; $HOURS_LOW / 20" | bc)
WEEKS_HIGH=$(echo "scale=1; $HOURS_HIGH / 20" | bc)
echo -e "${GREEN}Estimated Work Weeks:${NC} ${WEEKS_LOW}-${WEEKS_HIGH} weeks (@ 20hrs/week)"

echo ""
echo -e "${BLUE}========================================${NC}"
