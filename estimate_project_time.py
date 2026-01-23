#!/usr/bin/env python3
"""
Project Time Estimation Tool

Analyzes file modification times to estimate hands-on work time.
Assumes 8-hour workdays.

Usage:
    python3 estimate_project_time.py                              # Full analysis
    python3 estimate_project_time.py . 2025-11-07 2025-11-24      # Custom date range
"""

import os
import sys
from datetime import datetime
from collections import defaultdict
import argparse


# ANSI color codes
class Colors:
    BLUE = '\033[0;34m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    NC = '\033[0m'  # No Color


def format_header(text):
    """Format a header with color and separators."""
    sep = "=" * 78
    return f"{Colors.BLUE}{sep}\n{text}\n{sep}{Colors.NC}"


def collect_files_in_period(project_dir, start_ts, end_ts):
    """Collect all file modification times within the specified period during business hours (8am-6pm)."""
    files_data = []
    
    for root, dirs, filenames in os.walk(project_dir):
        # Skip hidden directories
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        
        for filename in filenames:
            path = os.path.join(root, filename)
            try:
                mtime = os.path.getmtime(path)
                if start_ts <= mtime <= end_ts:
                    # Only count modifications during business hours (8am-6pm)
                    mod_datetime = datetime.fromtimestamp(mtime)
                    if 8 <= mod_datetime.hour < 18:
                        files_data.append((mtime, path))
            except (OSError, FileNotFoundError):
                pass
    
    return sorted(files_data)


def count_file_extensions(files_data):
    """Count files by extension."""
    ext_count = {}
    for _, path in files_data:
        if '.' in path:
            ext = path.split('.')[-1]
            ext_count[ext] = ext_count.get(ext, 0) + 1
    return ext_count


def get_daily_breakdown(files_data):
    """Break down modifications by day."""
    days = defaultdict(int)
    for mtime, _ in files_data:
        day = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d')
        days[day] += 1
    return days


def format_bar_chart(count, total):
    """Create a simple bar chart for visualization."""
    if count >= 1000:
        return '█' * (count // 500)
    elif count >= 100:
        return '▓' * (count // 50)
    elif count >= 10:
        return '▒' * (count // 5)
    else:
        return '░'


def calculate_time_estimates(total_files):
    """Calculate work time estimates with different assumptions."""
    # Per-modification time estimates (in minutes)
    time_conservative = 0.5   # min
    time_moderate = 0.75       # min
    time_intensive = 1.25      # min
    
    # Calculate hours
    hours_conservative = (total_files * time_conservative) / 60
    hours_moderate = (total_files * time_moderate) / 60
    hours_intensive = (total_files * time_intensive) / 60
    
    # Calculate 8-hour workdays
    days_conservative = hours_conservative / 8
    days_moderate = hours_moderate / 8
    days_intensive = hours_intensive / 8
    
    # Calculate 5-day work weeks
    weeks_conservative = days_conservative / 5
    weeks_moderate = days_moderate / 5
    weeks_intensive = days_intensive / 5
    
    return {
        'conservative': {
            'hours': hours_conservative,
            'days': days_conservative,
            'weeks': weeks_conservative,
        },
        'moderate': {
            'hours': hours_moderate,
            'days': days_moderate,
            'weeks': weeks_moderate,
        },
        'intensive': {
            'hours': hours_intensive,
            'days': days_intensive,
            'weeks': weeks_intensive,
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description='Estimate project work time based on file modifications.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                           # Default analysis (Nov 7-24, 2025)
  %(prog)s --start-date 2025-11-07 --end-date 2025-11-24  # Custom date range
  %(prog)s -d /path/to/project                       # Custom directory
  %(prog)s -d /path/to/project -s 2025-01-01 -e 2025-12-31  # All options
        """
    )
    
    parser.add_argument('-d', '--directory', dest='project_dir', default='.',
                        help='Project directory (default: current directory)')
    parser.add_argument('-s', '--start-date', dest='start_date', default='2025-11-07',
                        help='Start date in YYYY-MM-DD format (default: 2025-11-07)')
    parser.add_argument('-e', '--end-date', dest='end_date', default='2025-11-24',
                        help='End date in YYYY-MM-DD format (default: 2025-11-24)')
    
    args = parser.parse_args()
    
    # Parse dates
    try:
        start_dt = datetime.strptime(args.start_date, '%Y-%m-%d')
        end_dt = datetime.strptime(args.end_date, '%Y-%m-%d').replace(
            hour=23, minute=59, second=59
        )
        start_ts = int(start_dt.timestamp())
        end_ts = int(end_dt.timestamp())
    except ValueError as e:
        print(f"Error: Invalid date format. Use YYYY-MM-DD\n{e}", file=sys.stderr)
        sys.exit(1)
    
    # Verify project directory exists
    if not os.path.isdir(args.project_dir):
        print(f"Error: Directory '{args.project_dir}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Collect files
    print(f"{Colors.BLUE}{'=' * 78}")
    print("PROJECT TIME ESTIMATION (File Modification Analysis)")
    print(f"{'=' * 78}{Colors.NC}\n")
    
    print(f"{Colors.GREEN}Project Directory:{Colors.NC} {args.project_dir}")
    print(f"{Colors.GREEN}Analysis Period:{Colors.NC}   {args.start_date} to {args.end_date}")
    print(f"{Colors.GREEN}Business Hours Filter:{Colors.NC} 8:00 AM - 6:00 PM\n")
    
    files_data = collect_files_in_period(args.project_dir, start_ts, end_ts)
    
    if not files_data:
        print("No files found in this period")
        sys.exit(0)
    
    # Calculate basic statistics
    earliest = files_data[0][0]
    latest = files_data[-1][0]
    total_files = len(files_data)
    span_hours = (latest - earliest) / 3600
    span_days = (latest - earliest) / 86400
    
    # Count active days
    active_dates = set()
    for mtime, _ in files_data:
        date = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d')
        active_dates.add(date)
    
    active_days = len(active_dates)
    
    # Print basic info
    print(f"Total file modifications:    {total_files:,}")
    print(f"Earliest modification:       {datetime.fromtimestamp(earliest).strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Latest modification:         {datetime.fromtimestamp(latest).strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Span:                        {span_days:.1f} calendar days ({span_hours:.0f} hours)")
    print(f"Active working days:         {active_days}\n")
    
    # Daily breakdown
    print(format_header("DAILY ACTIVITY"))
    days = get_daily_breakdown(files_data)
    
    for day in sorted(days.keys()):
        count = days[day]
        date_obj = datetime.strptime(day, '%Y-%m-%d')
        weekday = date_obj.strftime('%a')
        pct = (count / total_files) * 100
        bar = format_bar_chart(count, total_files)
        print(f"{day} ({weekday}): {count:7,d} mods ({pct:5.1f}%) {bar}")
    
    # File type breakdown
    print(f"\n{format_header('FILE TYPES')}\n")
    ext_count = count_file_extensions(files_data)
    
    for ext, count in sorted(ext_count.items(), key=lambda x: x[1], reverse=True)[:8]:
        pct = (count / total_files) * 100
        print(f"  .{ext:12s} {count:7,d} ({pct:5.1f}%)")
    
    # Time estimates
    print(f"\n{format_header('WORK TIME ESTIMATION (8-hour workdays)')}\n")
    
    estimates = calculate_time_estimates(total_files)
    
    print(f"{Colors.YELLOW}CONSERVATIVE (1.0 min per modification){Colors.NC}")
    print(f"  Total hours:      {estimates['conservative']['hours']:7.0f} hrs")
    print(f"  8-hour workdays:  {estimates['conservative']['days']:7.1f} days")
    print(f"  5-day work weeks: {estimates['conservative']['weeks']:7.1f} weeks\n")
    
    print(f"{Colors.YELLOW}MODERATE (1.5 min per modification){Colors.NC} {Colors.GREEN}← MOST LIKELY{Colors.NC}")
    print(f"  Total hours:      {estimates['moderate']['hours']:7.0f} hrs")
    print(f"  8-hour workdays:  {estimates['moderate']['days']:7.1f} days")
    print(f"  5-day work weeks: {estimates['moderate']['weeks']:7.1f} weeks\n")
    
    print(f"{Colors.YELLOW}INTENSIVE (2.5 min per modification, includes debugging){Colors.NC}")
    print(f"  Total hours:      {estimates['intensive']['hours']:7.0f} hrs")
    print(f"  8-hour workdays:  {estimates['intensive']['days']:7.1f} days")
    print(f"  5-day work weeks: {estimates['intensive']['weeks']:7.1f} weeks\n")
    
    # Per day averages
    print(f"{Colors.BLUE}{'=' * 78}")
    print(f"AVERAGE PER WORKING DAY ({active_days} days)")
    print(f"{'=' * 78}{Colors.NC}\n")
    
    avg_mods_per_day = total_files / active_days
    avg_hours_per_day = estimates['moderate']['hours'] / active_days
    
    print(f"Modifications per day:       {avg_mods_per_day:7,.0f}")
    print(f"Hours per day (moderate):    {min(avg_hours_per_day, 8):7.1f} hrs")
    
    print(f"  → {Colors.GREEN}8.0 hrs/day{Colors.NC} maximum (realistic daily capacity)")
    if avg_hours_per_day > 8:
        overtime_days = avg_hours_per_day / 8
        print(f"  → {Colors.YELLOW}{overtime_days:.1f} equivalent 8-hour days{Colors.NC} of work per calendar day")
    
    # Summary
    print(f"\n{Colors.BLUE}{'=' * 78}")
    print("SUMMARY")
    print(f"{'=' * 78}{Colors.NC}\n")
    
    hours_moderate = estimates['moderate']['hours']
    equivalent_workdays = hours_moderate / 8
    
    print(f"Total work hours: {Colors.GREEN}{hours_moderate:.0f} hours{Colors.NC}")
    print(f"Equivalent to: {Colors.GREEN}{equivalent_workdays:.1f} full 8-hour workdays{Colors.NC}")
    print(f"Compressed into: {Colors.GREEN}{active_days} actual calendar days{Colors.NC}\n")


if __name__ == '__main__':
    main()
