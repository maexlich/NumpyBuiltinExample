#include <boinc_api.h>
#include <app_ipc.h>

extern "C" int boinc_temporary_exit_wrapper(int delay, const char *reason, int is_notice){
	if (is_notice)
		return boinc_temporary_exit(delay,reason, true);
	else
		return boinc_temporary_exit(delay, reason, false);
}

extern "C" double boinc_get_ram(double fraction){
	APP_INIT_DATA init_data;
	HOST_INFO host_info;
	double ram_in_bytes = 0;
	boinc_get_init_data(init_data);
	host_info = init_data.host_info;
	ram_in_bytes = host_info.m_nbytes;
	return (ram_in_bytes/1073741824)*fraction;
}
