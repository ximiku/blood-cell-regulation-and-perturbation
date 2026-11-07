-- Official GraphPad AppleScript entry point; run only after user authorization.
-- Activates Prism and opens the task-specific native Prism batch script.
tell application "/Applications/Prism 10.app"
    activate
    open POSIX file "/path/to/user/Desktop/CourseData/26_fall/Statistical_analysis_of_genomics_data/HW2/solution/prism/migration/open_pca.pzc"
end tell
